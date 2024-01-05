// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphOutputModule.h"
#include "SMEffectDecoder.h"
#include "SMLeakDebugger.h"
#include "SMMorphOutput.h"
#include "SMMorphPlan.h"
#include "glib.h"
#include <assert.h>

#define CHANNEL_OP_COUNT 4

using namespace SpectMorph;

using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphOutputModule");

MorphOutputModule::MorphOutputModule(MorphPlanVoice* voice)
    : MorphOperatorModule(voice), decoder(this, morph_plan_voice->mix_freq()) {
    leak_debugger.add(this);
}

MorphOutputModule::~MorphOutputModule() {
    leak_debugger.del(this);
}

void MorphOutputModule::set_config(const MorphOperatorConfig* op_cfg) {
    cfg = dynamic_cast<const MorphOutput::Config*>(op_cfg);
    g_return_if_fail(cfg != nullptr);

    MorphOperatorModule* mod = morph_plan_voice->module(cfg->channel_ops[0]);
    LiveDecoderSource* source = mod ? mod->source() : nullptr;

    /* since the source is part of a module (and modules get newly created in
     * main thread and then replaced in audio thread), comparing the pointer to
     * the source in the LiveDecoder is enough to see if the source changed
     */
    decoder.set_config(cfg, source, morph_plan_voice->mix_freq());
}

void MorphOutputModule::process(const TimeInfoGenerator& time_info_gen_, RTMemoryArea& rt_memory_area, size_t n_samples,
                                float** values, const float* freq_in) {
    this->time_info_gen = &time_info_gen_;
    m_rt_memory_area = &rt_memory_area;

    decoder.process(rt_memory_area, n_samples, freq_in, values[0]);

    this->time_info_gen = nullptr;
    m_rt_memory_area = nullptr;
}

RTMemoryArea* MorphOutputModule::rt_memory_area() const {
    return m_rt_memory_area;
}

TimeInfo MorphOutputModule::compute_time_info() const {
    assert(time_info_gen);
    return time_info_gen->time_info(decoder.time_offset_ms());
}

void MorphOutputModule::prepareToPlay(float mix_freq) {
    decoder.prepareToPlay(mix_freq);
}

void MorphOutputModule::retrigger(const TimeInfo& time_info, int channel, float freq, int midi_velocity, bool onset) {
    decoder.retrigger(channel, freq, midi_velocity, onset);
}

void MorphOutputModule::release() {
    decoder.release();
}

bool MorphOutputModule::done() {
    // done means: the signal will be only zeros from here
    return decoder.done();
}
