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

MorphOutputModule::MorphOutputModule(MorphPlanVoice* voice) : MorphOperatorModule(voice) {
    out_ops.resize(CHANNEL_OP_COUNT);
    out_decoders.resize(CHANNEL_OP_COUNT);

    leak_debugger.add(this);
}

MorphOutputModule::~MorphOutputModule() {
    for (size_t ch = 0; ch < CHANNEL_OP_COUNT; ch++) {
        if (out_decoders[ch]) {
            delete out_decoders[ch];
            out_decoders[ch] = nullptr;
        }
    }
    leak_debugger.del(this);
}

void MorphOutputModule::set_config(const MorphOperatorConfig* op_cfg) {
    cfg = dynamic_cast<const MorphOutput::Config*>(op_cfg);
    g_return_if_fail(cfg != nullptr);

    clear_dependencies();
    for (size_t ch = 0; ch < CHANNEL_OP_COUNT; ch++) {
        EffectDecoder* dec = nullptr;

        MorphOperatorModule* mod = morph_plan_voice->module(cfg->channel_ops[ch]);

        if (mod == out_ops[ch]) // same source
        {
            dec = out_decoders[ch];
            // keep decoder as it is
        } else {
            if (out_decoders[ch])
                delete out_decoders[ch];
            if (mod) {
                dec = new EffectDecoder(mod->source());
            }
        }

        if (dec)
            dec->set_config(cfg);

        out_ops[ch] = mod;
        out_decoders[ch] = dec;

        add_dependency(mod);
    }
}

void MorphOutputModule::process(size_t n_samples, float** values, size_t n_ports, const float* freq_in) {
    g_return_if_fail(n_ports <= out_decoders.size());

    for (size_t port = 0; port < n_ports; port++) {
        if (values[port]) {
            if (out_decoders[port]) {
                out_decoders[port]->process(n_samples, freq_in, values[port]);
            } else {
                zero_float_block(n_samples, values[port]);
            }
        }
    }
}

void MorphOutputModule::prepareToPlay(float mix_freq) {
    for (size_t port = 0; port < CHANNEL_OP_COUNT; port++) {
        if (out_decoders[port]) {
            out_decoders[port]->prepareToPlay(mix_freq);
        }
    }
}

void MorphOutputModule::retrigger(int channel, float freq, int midi_velocity, bool onset) {
    for (size_t port = 0; port < CHANNEL_OP_COUNT; port++) {
        if (out_decoders[port]) {
            out_decoders[port]->retrigger(channel, freq, midi_velocity, onset);
        }
    }
}

void MorphOutputModule::release() {
    for (auto dec : out_decoders) {
        if (dec)
            dec->release();
    }
}
