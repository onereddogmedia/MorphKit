// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMModulationList.h"
#include "SMMorphPlan.h"

using namespace SpectMorph;

using std::string;
using std::vector;

ModulationList::ModulationList(ModulationData& data_, Property& property_) : data(data_), property(property_) {
    connect(property.op()->morph_plan()->signal_operator_removed, this, &ModulationList::on_operator_removed);
}

MorphOperator::ControlType ModulationList::main_control_type() const {
    return data.main_control_type;
}

MorphOperator* ModulationList::main_control_op() const {
    return data.main_control_op.get();
}

void ModulationList::set_main_control_type_and_op(MorphOperator::ControlType type, MorphOperator* op) {
    data.main_control_type = type;
    data.main_control_op.set(op);

    signal_main_control_changed();
    signal_modulation_changed();
}

size_t ModulationList::count() const {
    return data.entries.size();
}

void ModulationList::add_entry() {
    data.entries.emplace_back();

    signal_size_changed();
    signal_modulation_changed();
}

void ModulationList::update_entry(size_t index, ModulationData::Entry& new_entry) {
    data.entries[index] = new_entry;
    signal_modulation_changed();
}

const ModulationData::Entry& ModulationList::operator[](size_t index) const {
    return data.entries[index];
}

void ModulationList::remove_entry(size_t index) {
    g_return_if_fail(index >= 0 && index < data.entries.size());
    data.entries.erase(data.entries.begin() + (long)index);

    signal_size_changed();
    signal_modulation_changed();
}

void ModulationList::set_compat_type_and_op(const string& type, const string& op) {
    compat = true;

    compat_type_name = type;
    compat_op_name = op;
}

string ModulationList::event_name(const string& id, int index) {
    string s = property.identifier() + ".modulation." + id;

    if (index >= 0)
        s += string_printf("_%d", index);

    return s;
}

static bool starts_with(const string& key, const string& start) {
    return key.substr(0, start.size()) == start;
}

bool ModulationList::split_event_name(const string& name, const string& start, int& index) {
    string prefix = event_name(start) + "_";

    if (!starts_with(name, prefix))
        return false;

    index = atoi(name.substr(prefix.length()).c_str());
    return true;
}

void ModulationList::get_dependencies(vector<MorphOperator*>& deps) {
    if (data.main_control_type == MorphOperator::CONTROL_OP)
        deps.push_back(data.main_control_op.get());

    for (const auto& entry : data.entries)
        if (entry.control_type == MorphOperator::CONTROL_OP)
            deps.push_back(entry.control_op.get());
}

void ModulationList::on_operator_removed(MorphOperator* op) {
    // plan changed will be emitted automatically after remove, so we don't emit it here

    if (op == data.main_control_op.get()) {
        data.main_control_op.set(nullptr);
        if (data.main_control_type == MorphOperator::CONTROL_OP)
            data.main_control_type = MorphOperator::CONTROL_GUI;
    }
    uint index = 0;
    while (index < data.entries.size()) {
        if (op == data.entries[index].control_op.get()) {
            data.entries.erase(data.entries.begin() + index);
        } else {
            index++;
        }
    }
}
