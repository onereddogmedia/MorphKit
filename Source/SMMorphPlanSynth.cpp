// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphPlanSynth.h"
#include "SMLeakDebugger.h"
#include "SMMorphPlanVoice.h"

using namespace SpectMorph;

using std::map;
using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphPlanSynth");

MorphPlanSynth::MorphPlanSynth(size_t n_voices) : m_mix_freq(44100) {
    leak_debugger.add(this);

    m_wav_set_repo = new WavSetRepo();

    for (size_t i = 0; i < n_voices; i++)
        voices.push_back(new MorphPlanVoice(this));
}

MorphPlanSynth::~MorphPlanSynth() {
    leak_debugger.del(this);

    delete m_wav_set_repo;
    m_wav_set_repo = nullptr;

    for (size_t i = 0; i < voices.size(); i++) {
        delete voices[i];
        voices[i] = nullptr;
    }

    voices.clear();
}

MorphPlanVoice* MorphPlanSynth::voice(size_t i) const {
    g_return_val_if_fail(i < voices.size(), nullptr);

    return voices[i];
}

static vector<string> sorted_id_list(MorphPlanPtr plan) {
    vector<string> ids;

    if (plan) {
        const vector<MorphOperator*>& ops = plan->operators();
        for (vector<MorphOperator*>::const_iterator oi = ops.begin(); oi != ops.end(); oi++) {
            ids.push_back((*oi)->id());
        }
        sort(ids.begin(), ids.end());
    }
    return ids;
}

MorphPlanSynth::UpdateP MorphPlanSynth::prepare_update(MorphPlanPtr plan) /* main thread */
{
    UpdateP update = std::make_shared<Update>();

    for (auto o : plan->operators()) {
        MorphOperatorConfigP config(o->clone_config());
        update->new_configs.push_back(config);

        Update::Op op = {.ptr_id = o->ptr_id(), .type = o->type(), .config = config.get()};
        op.config->wav_set_repo = m_wav_set_repo;
        update->ops.push_back(op);
    }
    sort(update->ops.begin(), update->ops.end(),
         [](const Update::Op& a, const Update::Op& b) { return a.ptr_id < b.ptr_id; });

    vector<string> update_ids = sorted_id_list(plan);

    // TODO: cheap update doesn't make any sense
    update->cheap = (update_ids == m_last_update_ids) && (plan->id() == m_last_plan_id);
    m_last_update_ids = update_ids;
    m_last_plan_id = plan->id();

    return update;
}

void MorphPlanSynth::apply_update(MorphPlanSynth::UpdateP update) /* audio thread */
{
    /* life time for configs:
     *  - configs required for current update should be kept alive (m_active_configs)
     *  - configs no longer needed should be freed, but not in audio thread
     */
    update->old_configs = std::move(m_active_configs);
    m_active_configs = std::move(update->new_configs);

    if (update->cheap) {
        for (size_t i = 0; i < voices.size(); i++)
            voices[i]->cheap_update(update);
    } else {
        for (size_t i = 0; i < voices.size(); i++)
            voices[i]->full_update(update);
    }
}

void MorphPlanSynth::set_mix_freq(float mix_freq) {
    m_mix_freq = mix_freq;
    for (size_t i = 0; i < voices.size(); i++)
        voices[i]->set_mix_freq(mix_freq);
}

float MorphPlanSynth::mix_freq() const {
    return m_mix_freq;
}

bool MorphPlanSynth::have_output() const {
    if (voices.empty())
        return false;

    // all voices are the same: either each of them contains an output module, or none of them
    return voices[0]->output() != nullptr;
}
