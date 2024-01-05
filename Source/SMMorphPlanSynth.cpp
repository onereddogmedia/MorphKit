// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphPlanSynth.h"
#include "SMLeakDebugger.h"
#include "SMMorphOutputModule.h"
#include "SMMorphPlanVoice.h"

using namespace SpectMorph;

using std::map;
using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphPlanSynth");

MorphPlanSynth::MorphPlanSynth(size_t n_voices) {
    leak_debugger.add(this);

    for (size_t i = 0; i < n_voices; i++)
        voices.push_back(new MorphPlanVoice(m_mix_freq, this));
}

MorphPlanSynth::~MorphPlanSynth() {
    leak_debugger.del(this);

    for (size_t i = 0; i < voices.size(); i++)
        delete voices[i];

    voices.clear();
}

MorphPlanVoice* MorphPlanSynth::voice(size_t i) const {
    g_return_val_if_fail(i < voices.size(), nullptr);

    return voices[i];
}

static vector<string> sorted_id_list(const MorphPlan& plan) {
    vector<string> ids;

    for (auto op : plan.operators()) {
        ids.push_back(op->id());
    }
    sort(ids.begin(), ids.end());
    return ids;
}

static bool recursive_cycle_check(MorphOperator* start_op, int depth) {
    /* check if processing would fail due to cycles
     *
     * this check should avoid crashes in this situation, although no audio will be produced
     */
    if (depth > 500)
        return true;

    for (auto op : start_op->dependencies())
        if (op && recursive_cycle_check(op, depth + 1))
            return true;

    return false;
}

MorphPlanSynth::UpdateP MorphPlanSynth::prepare_update(const MorphPlan& plan) /* main thread */
{
    UpdateP update = std::make_shared<Update>();

    update->have_cycle = false;
    for (auto op : plan.operators()) {
        if (recursive_cycle_check(op, 0))
            update->have_cycle = true;
    }

    for (auto o : plan.operators()) {
        MorphOperatorConfig* config = o->clone_config();

        Update::Op op = {.ptr_id = o->ptr_id(), .type = o->type(), .config = config};
        update->ops.push_back(op);
        update->new_configs.emplace_back(config); // take ownership (unique_ptr)
    }
    sort(update->ops.begin(), update->ops.end(),
         [](const Update::Op& a, const Update::Op& b) { return a.ptr_id < b.ptr_id; });

    vector<string> update_ids = sorted_id_list(plan);

    update->cheap = (update_ids == m_last_update_ids) && (plan.id() == m_last_plan_id);
    m_last_update_ids = update_ids;
    m_last_plan_id = plan.id();
    if (!update->cheap) {
        update->voice_full_updates.resize(voices.size());
        update->new_shared_states.resize(update->ops.size());

        for (size_t voice = 0; voice < voices.size(); voice++) {
            for (size_t mod_index = 0; mod_index < update->ops.size(); mod_index++) {
                const auto& op = update->ops[mod_index];
                OpModule op_module;

                // avoid creating modules in audio thread by doing it here
                op_module.module.reset(MorphOperatorModule::create(op.type, voices[voice]));
                op_module.ptr_id = op.ptr_id;
                op_module.config = op.config;

                if (op_module.module) {
                    op_module.module->set_ptr_id(op.ptr_id);

                    /* setup one shared state structure for all voices */
                    if (voice == 0)
                        update->new_shared_states[mod_index].reset(op_module.module->create_shared_state());
                    if (update->new_shared_states[mod_index])
                        op_module.module->set_shared_state(update->new_shared_states[mod_index].get());

                    if (op.type == "SpectMorph::MorphOutput")
                        update->voice_full_updates[voice].output_module =
                            dynamic_cast<MorphOutputModule*>(op_module.module.get());

                    update->voice_full_updates[voice].new_modules.push_back(std::move(op_module));
                } else
                    g_warning("operator type %s lacks MorphOperatorModule\n", op.type.c_str());
            }
        }
    }

    return update;
}

void MorphPlanSynth::apply_update(MorphPlanSynth::UpdateP update) /* audio thread */
{
    /* life time for configs:
     *  - configs required for current update should be kept alive (m_active_configs)
     *  - configs no longer needed should be freed, but not in audio thread
     */
    m_active_configs.swap(update->new_configs);
    m_have_cycle = update->have_cycle;

    if (update->cheap) {
        for (size_t i = 0; i < voices.size(); i++)
            voices[i]->cheap_update(update);
    } else {
        voices_shared_states.swap(update->new_shared_states);

        for (size_t i = 0; i < voices.size(); i++)
            voices[i]->full_update(update->voice_full_updates[i]);
    }
}

void MorphPlanSynth::update_shared_state(const TimeInfo& time_info) {
    if (voices.empty())
        return;
    voices[0]->update_shared_state(time_info);
}

void MorphPlanSynth::set_mix_freq(float mix_freq) {
    m_mix_freq = mix_freq;
    for (size_t i = 0; i < voices.size(); i++)
        voices[i]->set_mix_freq(mix_freq);
}

float MorphPlanSynth::mix_freq() const {
    return m_mix_freq;
}

SpectMorph::Random* MorphPlanSynth::random_gen() {
    return &m_random_gen;
}

bool MorphPlanSynth::have_output() const {
    if (voices.empty())
        return false;

    // all voices are the same: either each of them contains an output module, or none of them
    return voices[0]->output() != nullptr;
}

bool MorphPlanSynth::have_cycle() const {
    return m_have_cycle;
}
