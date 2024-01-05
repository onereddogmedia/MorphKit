// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphSource.h"
#include "SMLeakDebugger.h"
#include "SMMorphPlan.h"

using namespace SpectMorph;

using std::string;

static LeakDebugger leak_debugger("SpectMorph::MorphSource");

MorphSource::MorphSource(MorphPlan* morph_plan) : MorphOperator(morph_plan) {
    leak_debugger.add(this);
}

MorphSource::~MorphSource() {
    leak_debugger.del(this);
}

void MorphSource::set_smset(const string& smset) {
    m_smset = smset;
    m_morph_plan->emit_plan_changed();
}

string MorphSource::smset() {
    return m_smset;
}

MorphOperatorConfig* MorphSource::clone_config() {
    Config* cfg = new Config(m_config);

    cfg->wav_set = m_morph_plan->wavSetRepo()->get(m_smset);

    return cfg;
}

const char* MorphSource::type() {
    return "SpectMorph::MorphSource";
}

int MorphSource::insert_order() {
    return 0;
}

MorphOperator::OutputType MorphSource::output_type() {
    return OUTPUT_AUDIO;
}
