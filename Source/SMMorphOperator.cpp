// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphOperator.h"
#include "SMMorphPlan.h"
#include "glib.h"

using namespace SpectMorph;

using std::map;
using std::string;
using std::vector;

MorphOperator::MorphOperator(MorphPlan* morph_plan) : m_morph_plan(morph_plan) {
}

MorphOperator::~MorphOperator() {
    // virtual destructor for proper subclass deletion
}

MorphPlan* MorphOperator::morph_plan() {
    return m_morph_plan;
}

string MorphOperator::name() {
    return m_name;
}

void MorphOperator::set_name(const string& name) {
    g_return_if_fail(can_rename(name));

    m_name = name;

    m_morph_plan->emit_plan_changed();
}

string MorphOperator::id() {
    return m_id;
}

void MorphOperator::set_id(const string& id) {
    m_id = id;
}

bool MorphOperator::can_rename(const string& name) {
    const vector<MorphOperator*>& ops = m_morph_plan->operators();

    if (name == "")
        return false;

    for (vector<MorphOperator*>::const_iterator oi = ops.begin(); oi != ops.end(); oi++) {
        MorphOperator* op = *oi;
        if (op != this && op->name() == name)
            return false;
    }
    return true;
}

string MorphOperator::type_name() {
    return string(type()).substr(string("SpectMorph::Morph").size());
}

vector<MorphOperator*> MorphOperator::dependencies() {
    return {}; /* default implementation -> no dependencies */
}

MorphOperatorConfig* MorphOperator::clone_config() {
    return nullptr; // FIXME: remove default impl
}

MorphOperatorConfig::~MorphOperatorConfig() {
}

#include "SMMorphGrid.h"
#include "SMMorphOutput.h"

MorphOperator* MorphOperator::create(const string& type, MorphPlan* plan) {
    g_return_val_if_fail(plan != nullptr, nullptr);

    if (type == "SpectMorph::MorphGrid")
        return new MorphGrid(plan);
    if (type == "SpectMorph::MorphOutput")
        return new MorphOutput(plan);

    return nullptr;
}
