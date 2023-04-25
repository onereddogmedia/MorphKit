// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphPlanVoice.h"
#include "SMEffectDecoder.h"
#include "SMLeakDebugger.h"
#include "SMMorphOutput.h"
#include "SMMorphOutputModule.h"
#include <assert.h>
#include <map>

using namespace SpectMorph;
using std::map;
using std::max;
using std::sort;
using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphPlanVoice");

MorphPlanVoice::MorphPlanVoice(MorphPlanSynth* synth)
    : m_control_input(MorphPlan::N_CONTROL_INPUTS), m_control_output(MorphPlan::N_CONTROL_INPUTS), m_output(nullptr),
      m_mix_freq(44100), m_morph_plan_synth(synth) {
    leak_debugger.add(this);
}

void MorphPlanVoice::configure_modules() {
    for (size_t i = 0; i < modules.size(); i++)
        modules[i].module->set_config(modules[i].config);
}

void MorphPlanVoice::create_modules(MorphPlanSynth::UpdateP update) {
    for (auto& op : update->ops) {
        MorphOperatorModule* module = MorphOperatorModule::create(op.type, this);

        if (!module) {
            g_warning("operator type %s lacks MorphOperatorModule\n", op.type.c_str());
        } else {
            module->set_ptr_id(op.ptr_id);

            OpModule op_module;

            op_module.module = module;
            op_module.ptr_id = op.ptr_id;
            op_module.config = op.config;

            modules.push_back(op_module);

            if (op.type == "SpectMorph::MorphOutput")
                m_output = dynamic_cast<MorphOutputModule*>(module);
        }
    }
}

void MorphPlanVoice::clear_modules() {
    for (size_t i = 0; i < modules.size(); i++) {
        assert(modules[i].module != nullptr);
        delete modules[i].module;
        modules[i].module = nullptr;
    }
    modules.clear();

    m_output = nullptr;
}

MorphPlanVoice::~MorphPlanVoice() {
    clear_modules();
    leak_debugger.del(this);
}

MorphOutputModule* MorphPlanVoice::output() {
    return m_output;
}

MorphOperatorModule* MorphPlanVoice::module(const MorphOperatorPtr& ptr) {
    MorphOperator::PtrID ptr_id = ptr.ptr_id();

    for (size_t i = 0; i < modules.size(); i++)
        if (modules[i].ptr_id == ptr_id)
            return modules[i].module;

    return nullptr;
}

void MorphPlanVoice::full_update(MorphPlanSynth::UpdateP update) {
    /* This will loose the original state information which means the audio
     * will not transition smoothely. However, this should only occur for plan
     * changes, not parameter updates.
     */
    clear_modules();
    create_modules(update);
    configure_modules();
}

void MorphPlanVoice::cheap_update(MorphPlanSynth::UpdateP update) {
    g_return_if_fail(update->ops.size() == modules.size());

    // exchange old operators with new operators
    for (size_t i = 0; i < modules.size(); i++) {
        assert(modules[i].ptr_id == update->ops[i].ptr_id);
        modules[i].config = update->ops[i].config;
        assert(modules[i].config);
    }

    // reconfigure modules
    configure_modules();
}

double MorphPlanVoice::control_input(double value, MorphOperator::ControlType ctype) {
    switch (ctype) {
        case MorphOperator::CONTROL_GUI:
            return value;
        case MorphOperator::CONTROL_SIGNAL_1:
            return m_control_input[0];
        case MorphOperator::CONTROL_SIGNAL_2:
            return m_control_input[1];
        case MorphOperator::CONTROL_SIGNAL_3:
            return m_control_input[2];
        case MorphOperator::CONTROL_SIGNAL_4:
            return m_control_input[3];
        case MorphOperator::CONTROL_SIGNAL_5:
            return m_control_input[4];
        case MorphOperator::CONTROL_SIGNAL_6:
            return m_control_input[5];
        case MorphOperator::CONTROL_OP:
        default:
            g_assert_not_reached();
            return 0;
    }
}

void MorphPlanVoice::set_control_input(unsigned int i, double value) {
    assert(i >= 0 && i < MorphPlan::N_CONTROL_INPUTS);

    m_control_input[i] = value;
}

double MorphPlanVoice::control_output(unsigned int i) {
    assert(i >= 0 && i < MorphPlan::N_CONTROL_INPUTS);
    return m_control_output[i];
}

void MorphPlanVoice::set_control_output(unsigned int i, double value) {
    assert(i >= 0 && i < MorphPlan::N_CONTROL_INPUTS);

    m_control_output[i] = value;
}

void MorphPlanVoice::set_mix_freq(float freq) {
    m_mix_freq = freq;
    for (size_t i = 0; i < modules.size(); i++) {
        modules[i].module->prepareToPlay(freq);
    }
    if (m_output)
        m_output->prepareToPlay(freq);
}

float MorphPlanVoice::mix_freq() const {
    return m_mix_freq;
}

MorphPlanSynth* MorphPlanVoice::morph_plan_synth() const {
    return m_morph_plan_synth;
}
