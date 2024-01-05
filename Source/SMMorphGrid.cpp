// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphGrid.h"
#include "SMLeakDebugger.h"
#include "SMMorphPlan.h"
#include "SMUtils.h"
#include "SMWavSet.h"
#include <assert.h>

using namespace SpectMorph;

using std::make_pair;
using std::map;
using std::pair;
using std::string;
using std::u32string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphGrid");

MorphGrid::MorphGrid(MorphPlan* morph_plan) : MorphOperator(morph_plan) {
    leak_debugger.add(this);

    connect(morph_plan->signal_operator_removed, this, &MorphGrid::on_operator_removed);

    m_config.width = 2;
    m_config.height = 1;

    update_size();
}

MorphGrid::~MorphGrid() {
    leak_debugger.del(this);
}

const char* MorphGrid::type() {
    return "SpectMorph::MorphGrid";
}

int MorphGrid::insert_order() {
    return 700;
}

MorphOperator::OutputType MorphGrid::output_type() {
    return OUTPUT_AUDIO;
}

void MorphGrid::set_width(int width) {
    m_config.width = width;
    update_size();

    m_morph_plan->emit_plan_changed();
}

int MorphGrid::width() {
    return m_config.width;
}

void MorphGrid::set_height(int height) {
    m_config.height = height;
    update_size();

    m_morph_plan->emit_plan_changed();
}

int MorphGrid::height() {
    return m_config.height;
}

void MorphGrid::update_size() {
    m_config.input_node.resize((size_t)m_config.width);
    for (int i = 0; i < m_config.width; i++)
        m_config.input_node[(size_t)i].resize((size_t)m_config.height);
}

MorphGridNode MorphGrid::input_node(int x, int y) {
    g_return_val_if_fail(x >= 0 && x < m_config.width, MorphGridNode());
    g_return_val_if_fail(y >= 0 && y < m_config.height, MorphGridNode());

    return m_config.input_node[(size_t)x][(size_t)y];
}

void MorphGrid::set_input_node(int x, int y, const MorphGridNode& node) {
    g_return_if_fail(x >= 0 && x < m_config.width);
    g_return_if_fail(y >= 0 && y < m_config.height);
    g_return_if_fail(node.smset == "" || !node.op); // should not set both

    m_config.input_node[(size_t)x][(size_t)y] = node;
    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_x_morphing(double new_morphing) {
    m_config.x_morphing = new_morphing;
}

void MorphGrid::set_y_morphing(double new_morphing) {
    m_config.y_morphing = new_morphing;
}

void MorphGrid::on_operator_removed(MorphOperator* op) {
    // plan changed will be emitted automatically after remove, so we don't emit it here

    for (size_t x = 0; x < (size_t)m_config.width; x++) {
        for (size_t y = 0; y < (size_t)m_config.height; y++) {
            if (m_config.input_node[x][y].op.get() == op) {
                assert(m_config.input_node[x][y].smset.empty());

                m_config.input_node[x][y].op.set(nullptr);
            }
        }
    }
}

vector<MorphOperator*> MorphGrid::dependencies() {
    std::vector<MorphOperator*> deps;

    for (size_t x = 0; x < (size_t)m_config.width; x++) {
        for (size_t y = 0; y < (size_t)m_config.height; y++) {
            deps.push_back(m_config.input_node[x][y].op.get());
        }
    }
    return deps;
}

MorphOperatorConfig* MorphGrid::clone_config() {
    Config* cfg = new Config(m_config);

    for (size_t x = 0; x < (size_t)m_config.width; x++) {
        for (size_t y = 0; y < (size_t)m_config.height; y++) {
            string smset = cfg->input_node[x][y].smset;
            if (smset != "") {
                cfg->input_node[x][y].wav_set = m_morph_plan->wavSetRepo()->get(smset);
            } else {
                cfg->input_node[x][y].wav_set = nullptr;
            }
        }
    }
    return cfg;
}

MorphGridNode::MorphGridNode() {
}
