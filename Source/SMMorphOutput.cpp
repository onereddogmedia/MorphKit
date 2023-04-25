// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphOutput.h"
#include "SMLeakDebugger.h"
#include "SMMorphPlan.h"

#include <assert.h>

#define CHANNEL_OP_COUNT 4

using namespace SpectMorph;

using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphOutput");

MorphOutput::MorphOutput(MorphPlan* morph_plan) : MorphOperator(morph_plan) {
    connect(morph_plan->signal_operator_removed, this, &MorphOutput::on_operator_removed);

    m_config.channel_ops.resize(CHANNEL_OP_COUNT);

    leak_debugger.add(this);
}

MorphOutput::~MorphOutput() {
    leak_debugger.del(this);
}

const char* MorphOutput::type() {
    return "SpectMorph::MorphOutput";
}

int MorphOutput::insert_order() {
    return 1000;
}

MorphOperator::OutputType MorphOutput::output_type() {
    return OUTPUT_NONE;
}

void MorphOutput::set_channel_op(int ch, MorphOperator* op) {
    assert(ch >= 0 && ch < CHANNEL_OP_COUNT);

    m_config.channel_ops[(size_t)ch].set(op);
    m_morph_plan->emit_plan_changed();
}

MorphOperator* MorphOutput::channel_op(int ch) {
    assert(ch >= 0 && ch < CHANNEL_OP_COUNT);

    return m_config.channel_ops[(size_t)ch].get();
}

void MorphOutput::on_operator_removed(MorphOperator* op) {
    for (size_t ch = 0; ch < m_config.channel_ops.size(); ch++) {
        if (m_config.channel_ops[ch].get() == op)
            m_config.channel_ops[ch].set(nullptr);
    }
}

vector<MorphOperator*> MorphOutput::dependencies() {
    vector<MorphOperator*> deps;

    for (auto& ptr : m_config.channel_ops)
        deps.push_back(ptr.get());

    return deps;
}

MorphOperatorConfig* MorphOutput::clone_config() {
    Config* cfg = new Config(m_config);
    return cfg;
}
