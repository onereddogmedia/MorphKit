// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphOperatorModule.h"

#include "SMLeakDebugger.h"
#include "SMMorphGridModule.h"
#include "SMMorphOutputModule.h"
#include "SMMorphPlanSynth.h"

using namespace SpectMorph;

using std::string;
using std::vector;

MorphOperatorModule::MorphOperatorModule(MorphPlanVoice* voice) : morph_plan_voice(voice) {
}

MorphOperatorModule::~MorphOperatorModule() {
    // virtual destructor to allow subclass deletion
}

void MorphOperatorModule::prepareToPlay(float mix_freq) {
    if (source())
        source()->prepareToPlay(mix_freq);
}

LiveDecoderSource* MorphOperatorModule::source() {
    return nullptr;
}

float MorphOperatorModule::value() {
    return 0;
}

void MorphOperatorModule::clear_dependencies() {
    m_dependencies.clear();
}

void MorphOperatorModule::add_dependency(MorphOperatorModule* dep_mod) {
    if (dep_mod)
        m_dependencies.push_back(dep_mod);
}

const vector<MorphOperatorModule*>& MorphOperatorModule::dependencies() const {
    return m_dependencies;
}

int& MorphOperatorModule::update_value_tag() {
    return m_update_value_tag;
}

void MorphOperatorModule::set_ptr_id(MorphOperator::PtrID ptr_id) {
    m_ptr_id = ptr_id;
}

MorphOperatorModule* MorphOperatorModule::create(const std::string& type, MorphPlanVoice* voice) {
    if (type == "SpectMorph::MorphGrid")
        return new MorphGridModule(voice);
    if (type == "SpectMorph::MorphOutput")
        return new MorphOutputModule(voice);

    return nullptr;
}
