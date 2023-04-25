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

    m_config.width = 2;
    m_config.height = 1;
    m_config.x_morphing = 0;
    m_config.y_morphing = 0;
    m_config.x_control_type = CONTROL_GUI;
    m_config.y_control_type = CONTROL_GUI;
    m_config.node_a_db_control_type = CONTROL_GUI;
    m_config.node_b_db_control_type = CONTROL_GUI;
    m_config.node_c_db_control_type = CONTROL_GUI;
    m_config.node_d_db_control_type = CONTROL_GUI;

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

void MorphGrid::set_width(unsigned int width) {
    m_config.width = width;
    update_size();

    m_morph_plan->emit_plan_changed();
}

unsigned int MorphGrid::width() {
    return m_config.width;
}

void MorphGrid::set_height(unsigned int height) {
    m_config.height = height;
    update_size();

    m_morph_plan->emit_plan_changed();
}

unsigned int MorphGrid::height() {
    return m_config.height;
}

void MorphGrid::update_size() {
    m_config.input_node.resize(m_config.width);
    for (unsigned int i = 0; i < m_config.width; i++)
        m_config.input_node[i].resize(m_config.height);
}

MorphGridNode MorphGrid::input_node(unsigned int x, unsigned int y) {
    g_return_val_if_fail(x < m_config.width, MorphGridNode());
    g_return_val_if_fail(y < m_config.height, MorphGridNode());

    return m_config.input_node[x][y];
}

void MorphGrid::set_input_node(unsigned int x, unsigned int y, const MorphGridNode& node) {
    g_return_if_fail(x < m_config.width);
    g_return_if_fail(y < m_config.height);

    m_config.input_node[x][y] = node;
    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_x_morphing(double new_morphing) {
    m_config.x_morphing = new_morphing;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_x_control_type(MorphGrid::ControlType control_type) {
    m_config.x_control_type = control_type;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_y_morphing(double new_morphing) {
    m_config.y_morphing = new_morphing;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_y_control_type(MorphGrid::ControlType control_type) {
    m_config.y_control_type = control_type;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_node_a_db_control_type(MorphGrid::ControlType control_type) {
    m_config.node_a_db_control_type = control_type;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_node_b_db_control_type(MorphGrid::ControlType control_type) {
    m_config.node_b_db_control_type = control_type;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_node_c_db_control_type(MorphGrid::ControlType control_type) {
    m_config.node_c_db_control_type = control_type;

    m_morph_plan->emit_plan_changed();
}

void MorphGrid::set_node_d_db_control_type(MorphGrid::ControlType control_type) {
    m_config.node_d_db_control_type = control_type;

    m_morph_plan->emit_plan_changed();
}

MorphOperatorConfig* MorphGrid::clone_config() {
    Config* cfg = new Config(m_config);

    for (unsigned int x = 0; x < m_config.width; x++) {
        for (unsigned int y = 0; y < m_config.height; y++) {
            string smset = cfg->input_node[x][y].smset;
            if (smset != "") {
                cfg->input_node[x][y].path = "/" + smset;
            } else {
                cfg->input_node[x][y].path = "";
            }
        }
    }
    return cfg;
}

MorphGridNode::MorphGridNode() {
}
