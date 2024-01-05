// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "MorphKit.h"
#include "SMAudioTool.h"
#include "SMEncoder.h"
#include "SMFft.h"
#include "SMLeakDebugger.h"
#include "SMWavData.h"

using namespace SpectMorph;

static LeakDebugger leak_debugger("SpectMorph::MorphKit");

MorphKit::MorphKit(const std::string& path, const std::string& icloud) : m_time_info_gen(48000) {
    SpectMorph::sm_set_pkg_data_dir(path);
    SpectMorph::sm_set_icloud_data_dir(icloud);
    int_sincos_init();
    sm_math_init();
    project.reset(new SpectMorph::Project());

    project->set_mix_freq(48000);
    makeMorphPlan();
}

MorphKit::~MorphKit() {
}

void MorphKit::setSampleRate(const double newRate) {
    project->set_mix_freq(newRate);
}

void MorphKit::setModel(const uint model, const std::string& name) {
    for (auto op : project->morph_plan()->operators()) {
        std::string type = op->type();
        if (type == "SpectMorph::MorphGrid") {
            int x = 0;
            int y = 0;
            if (model == 1) { // A
                x = 0;
                y = 1;
            } else if (model == 2) { // D
                x = 0;
                y = 0;
            } else if (model == 3) { // B
                x = 1;
                y = 1;
            } else if (model == 4) { // C
                x = 1;
                y = 0;
            }
            SpectMorph::MorphGrid* morph_grid = dynamic_cast<SpectMorph::MorphGrid*>(op);
            SpectMorph::MorphGridNode node = morph_grid->input_node(x, y);
            node.smset = name;
            morph_grid->set_input_node(x, y, node);
        }
    }

    project->morph_plan()->signal_plan_changed();
    project->try_update_synth();
}

void MorphKit::setLinear(const float value) {
    for (auto op : project->morph_plan()->operators()) {
        std::string type = op->type();
        if (type == "SpectMorph::MorphGrid") {
            SpectMorph::MorphGrid* morph_grid = dynamic_cast<SpectMorph::MorphGrid*>(op);
            int width;
            int height = 2;
            if (value < 2) {
                width = 2;
            } else if (approximatelyEqual(value, 2.0f)) {
                width = 1;
            } else {
                width = 1;
                height = 1;
            }
            morph_grid->set_width(width);
            morph_grid->set_height(height);

            SpectMorph::MorphGridNode node = morph_grid->input_node(0, 1);
            morph_grid->set_input_node(0, 1, node);
            node = morph_grid->input_node(0, 0);
            morph_grid->set_input_node(0, 0, node);
            node = morph_grid->input_node(1, 1);
            morph_grid->set_input_node(1, 1, node);
            node = morph_grid->input_node(1, 0);
            morph_grid->set_input_node(1, 0, node);
        }
    }
}

void MorphKit::reloadWavSet() {
    project->morph_plan()->reloadWavSet();
    project->try_update_synth();
}

void MorphKit::makeMorphPlan() {
    project->set_state_changed_notify(false);
    auto plan = project->morph_plan();

    auto grid = dynamic_cast<SpectMorph::MorphGrid*>(SpectMorph::MorphOperator::create("SpectMorph::MorphGrid", plan));
    grid->set_id(plan->generate_id());
    grid->set_name("Grid #1");
    grid->set_width(2);
    grid->set_height(2);
    grid->set_x_morphing(0);
    grid->set_y_morphing(0);

    SpectMorph::MorphGridNode node1;
    node1.smset = "/Synth/Sine.smset";
    SpectMorph::MorphGridNode node2;
    node2.smset = "/Synth/Sine.smset";
    SpectMorph::MorphGridNode node3;
    node3.smset = "/Synth/Sine.smset";
    SpectMorph::MorphGridNode node4;
    node4.smset = "/Synth/Sine.smset";

    /*
     *  A ---- B
     *  |      |
     *  |      |
     *  C ---- D
     */
    grid->set_input_node(0, 0, node1); // D
    grid->set_input_node(1, 0, node2); // C
    grid->set_input_node(0, 1, node3); // A
    grid->set_input_node(1, 1, node4); // B

    plan->add_operator(grid, SpectMorph::MorphPlan::ADD_POS_END);

    auto out =
        dynamic_cast<SpectMorph::MorphOutput*>(SpectMorph::MorphOperator::create("SpectMorph::MorphOutput", plan));
    out->set_channel_op(0, grid);
    plan->add_operator(out, SpectMorph::MorphPlan::ADD_POS_END);

    project->set_state_changed_notify(true);
    project->try_update_synth();
}

void MorphKitVoice::startNote(const float note_frequency, const int8_t midi_velocity, bool onset) {
    TimeInfo time_info = m_time_info_gen->time_info(0);
    if (mp_voice) {
        SpectMorph::MorphOutputModule* output_module = mp_voice->output();
        if (output_module) {
            output_module->retrigger(time_info, 0, note_frequency, midi_velocity, onset);
            mp_voice->set_control_input(0, 0);
            mp_voice->set_control_input(1, 0);
        }
    }
}

void MorphKitVoice::setControl(const double x, const double y, const double gain_a, const double gain_b,
                               const double gain_c, const double gain_d) {
    if (mp_voice) {
        mp_voice->set_control_input(0, x);
        mp_voice->set_control_input(1, y);
        mp_voice->set_control_input(2, gain_a);
        mp_voice->set_control_input(3, gain_b);
        mp_voice->set_control_input(4, gain_c);
        mp_voice->set_control_input(5, gain_d);
    }
}

void MorphKitVoice::prepareToPlay(const size_t voice_id) {
    if (morph_plan) {
        mp_voice = morph_plan->project()->midi_synth()->plan_synth()->voice(voice_id);
    }
}

template <typename sampleType>
void MorphKitVoice::render(const double frequency, sampleType* out, const size_t numSamples) {
    float morph[kMaxBlockSize];

    if (mp_voice) {
        SpectMorph::MorphOutputModule* output_module = mp_voice->output();
        if (output_module) {
            float* values[1] = {morph};
            float frequencies[kMaxBlockSize];
            for (size_t sample = 0; sample < numSamples; ++sample) {
                frequencies[sample] = (float)frequency;
            }
            output_module->process(*m_time_info_gen, *rt_memory_area, numSamples, values, frequencies);
            for (size_t sample = 0; sample < numSamples; ++sample) {
                out[sample] += morph[sample];
            }
        }
    }
}

// explicit instantiation for supported float types:
template void MorphKitVoice::render(const double frequency, float* out, const size_t numSamples);
template void MorphKitVoice::render(const double frequency, double* out, const size_t numSamples);

void MorphKitVoice::getRMS(float& rms_a, float& rms_b, float& rms_c, float& rms_d) {
    if (mp_voice) {
        rms_a = (float)mp_voice->control_output(0);
        rms_b = (float)mp_voice->control_output(1);
        rms_c = (float)mp_voice->control_output(2);
        rms_d = (float)mp_voice->control_output(3);
    }
}

void MorphKitEncoder::encodeToFile(std::vector<float>& samples, int channels, float sampleRate, int bitsPerSample,
                                   int midi_note, float fundamental_freq, int start, int end, Audio::LoopType loop_type,
                                   int loop_start, int loop_end, bool auto_vol, float& norm_db, bool auto_pitch,
                                   const std::string& sample_name, const std::string& path,
                                   std::function<bool()> kill_function) {
    WavData wav_data(samples, channels, sampleRate, bitsPerSample);

    const int optimization_level = 1;
    const bool track_sines = true;

    SpectMorph::EncoderParams enc_params;
    enc_params.setup_params(wav_data, fundamental_freq);
    std::vector<float> window(enc_params.block_size);

    std::string window_type;
    if (!enc_params.get_param("window", window_type))
        window_type = "hann";

    for (uint i = 0; i < window.size(); i++) {
        const size_t frame_size = enc_params.frame_size;

        if (i < frame_size) {
            if (window_type == "hann") {
                // hann
                window[i] = (float)window_cos(2.0 * i / (frame_size - 1) - 1.0);
            } else if (window_type == "hamming") {
                // probably never a good idea, since the sidelobes of the spectrum
                // do not roll off fast (as with the hann window)
                window[i] = (float)window_hamming(2.0 * i / (frame_size - 1) - 1.0);
            } else if (window_type == "blackman") {
                window[i] = (float)window_blackman(2.0 * i / (frame_size - 1) - 1.0);
            } else {
                printf("unsupported window type in config.\n");
            }
        } else {
            window[i] = 0;
        }
    }
    enc_params.window = window;
    enc_params.set_kill_function(kill_function);
    Encoder encoder(enc_params);
    bool result = encoder.encode(wav_data, 0, optimization_level, track_sines);
    if (!result) {
        return;
    }

    // strip models
    std::vector<SpectMorph::EncoderBlock>& audio_blocks = encoder.audio_blocks;
    for (size_t i = 0; i < audio_blocks.size(); i++) {
        audio_blocks[i].debug_samples.clear();
    }

    encoder.set_loop(loop_type, loop_start, loop_end);

    Audio* audio = encoder.save_as_audio();
    const int last_frame = audio->contents.size() ? ((const int)audio->contents.size() - 1) : 0;
    const double zero_values_ms = audio->zero_values_at_start / audio->mix_freq * 1000.0;
    double start_ms = start * 1000.0 / sampleRate;
    double end_ms = end * 1000.0 / sampleRate;
    double loop_start_ms = loop_start * 1000.0 / sampleRate;
    double loop_end_ms = loop_end * 1000.0 / sampleRate;
    const int x_start = sm_bound<int>(0, (int)lrint((zero_values_ms + start_ms) / audio->frame_step_ms), last_frame);
    const int x_end = sm_bound<int>(0, (int)lrint((zero_values_ms + end_ms) / audio->frame_step_ms), last_frame);
    const int xloop_start =
        sm_bound<int>(0, (int)lrint((zero_values_ms + loop_start_ms) / audio->frame_step_ms), last_frame);
    const int xloop_end =
        sm_bound<int>(0, (int)lrint((zero_values_ms + loop_end_ms) / audio->frame_step_ms), last_frame);
    audio->start = x_start;
    audio->end = x_end;
    audio->loop_start = xloop_start;
    audio->loop_end = xloop_end;
    audio->original_start = start;
    audio->original_end = end;
    audio->original_loop_start = loop_start;
    audio->original_loop_end = loop_end;

    if (auto_pitch) {
        applyAutoTune(*audio);
    }
    if (auto_vol) {
        applyAutoVolume(*audio);
        norm_db = audio->original_samples_norm_db;
    }

    WavSet wav_set;
    wav_set.name = sample_name;
    wav_set.short_name = wav_set.name;
    WavSetWave new_wave;
    new_wave.midi_note = midi_note;
    new_wave.channel = 0;
    new_wave.velocity_range_min = 0;
    new_wave.velocity_range_max = 127;
    new_wave.audio = audio;
    wav_set.waves.push_back(new_wave);
    wav_set.save(path, true);
}

void MorphKitEncoder::applyAutoVolume(Audio& audio) {
    double energy = AudioTool::compute_energy(audio);
    AudioTool::normalize_energy(energy, audio);
}

void MorphKitEncoder::applyAutoTune(Audio& audio) {
    // TODO: config autotune
    struct AutoTune {
        enum { SIMPLE, ALL_FRAMES, SMOOTH } method = SIMPLE;
        bool enabled = false;
        int partials = 1;   // used by: all_frames, smooth
        double time = 100;  // used_by: smooth
        double amount = 25; // used by: smooth
    };
    AutoTune auto_tune;
    auto_tune.method = AutoTune::SIMPLE;
    if (auto_tune.method == AutoTune::SIMPLE) {
        double tune_factor;
        if (AudioTool::get_auto_tune_factor(audio, tune_factor))
            AudioTool::apply_auto_tune_factor(audio, tune_factor);
    } else if (auto_tune.method == AutoTune::ALL_FRAMES) {
        for (auto& block : audio.contents) {
            const double est_freq = block.estimate_fundamental(auto_tune.partials);
            const double tune_factor = 1.0 / est_freq;

            for (size_t p = 0; p < block.freqs.size(); p++) {
                const double freq = block.freqs_f(p) * tune_factor;
                block.freqs[p] = SpectMorph::sm_freq2ifreq(freq);
            }
        }
    } else if (auto_tune.method == AutoTune::SMOOTH) {
        AudioTool::auto_tune_smooth(audio, auto_tune.partials, auto_tune.time, auto_tune.amount);
    }
}
