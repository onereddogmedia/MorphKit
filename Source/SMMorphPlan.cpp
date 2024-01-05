// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphPlan.h"
#include "SMAudio.h"
#include "SMInFile.h"
#include "SMLeakDebugger.h"
#include "SMMMapIn.h"
#include "SMMemOut.h"
#include "SMMorphOutput.h"
#include "SMOutFile.h"
#include "SMUtils.h"

#include <assert.h>
#include <map>
#include <random>

using namespace SpectMorph;

using std::map;
using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphPlan");

MorphPlan::MorphPlan(Project& project) : m_project(&project) {
    in_restore = false;
    m_id = generate_id();

    m_wav_set_repo = new WavSetRepo();

    leak_debugger.add(this);
}

/**
 * Clears MorphPlan, deletes all operators and resets all variables to return
 * to a state that is completely empty (like after creation).
 */
void MorphPlan::clear() {
    for (vector<MorphOperator*>::iterator oi = m_operators.begin(); oi != m_operators.end(); oi++)
        delete (*oi);
    m_operators.clear();

    /*
     * generate a new MorphPlan id; the id is used by MorphPlanSynth to detect the case
     * that after a reload of the same .smplan, all operators have the same id, but still
     * a full update of the modules is necessary
     */
    m_id = generate_id();
}

MorphPlan::~MorphPlan() {
    assert(!in_restore);

    clear();

    delete m_wav_set_repo;
    m_wav_set_repo = nullptr;

    leak_debugger.del(this);
}

void MorphPlan::add_operator(MorphOperator* op, AddPos add_pos, const string& load_name, const string& load_id,
                             bool load_folded) {
    if (load_name == "") {
        // generate uniq name
        string name;
        bool uniq;
        int i = 0;

        do {
            i++;
            uniq = true;
            name = string_printf("%s #%d", op->type_name().c_str(), i);

            for (vector<MorphOperator*>::iterator oi = m_operators.begin(); oi != m_operators.end(); oi++) {
                if ((*oi)->name() == name)
                    uniq = false;
            }
        } while (!uniq);
        op->set_name(name);
    } else {
        op->set_name(load_name);
    }
    if (load_id == "") {
        op->set_id(generate_id());
    } else {
        op->set_id(load_id);
    }
    op->set_folded(load_folded);

    if (add_pos == ADD_POS_AUTO) {
        size_t pos = 0;

        for (size_t i = 0; i < m_operators.size(); i++) {
            if (m_operators[i]->insert_order() <= op->insert_order())
                pos = i + 1;
        }
        m_operators.insert(m_operators.begin() + (long)pos, op);
    } else {
        m_operators.push_back(op);
    }

    signal_operator_added(op);
    signal_need_view_rebuild();
    emit_plan_changed();
}

void MorphPlan::reloadWavSet() {
    if (m_wav_set_repo) {
        m_wav_set_repo->reload();
    }
    emit_plan_changed();
}

/**
 * Get MorphPlan operators.
 *
 * \returns a read-only reference to the vector containing the operators.
 */
const vector<MorphOperator*>& MorphPlan::operators() const {
    return m_operators;
}

void MorphPlan::remove(MorphOperator* op) {
    signal_need_view_rebuild();
    signal_operator_removed(op);

    // accessing operator contents after remove was called is an error
    delete op;

    vector<MorphOperator*>::iterator oi = m_operators.begin();
    while (oi != m_operators.end()) {
        if (*oi == op)
            oi = m_operators.erase(oi);
        else
            oi++;
    }

    emit_plan_changed();
}

void MorphPlan::move(MorphOperator* op, MorphOperator* op_next) {
    signal_need_view_rebuild();

    vector<MorphOperator*> new_operators;

    // change op position so that op is before op_next
    for (vector<MorphOperator*>::iterator oi = m_operators.begin(); oi != m_operators.end(); oi++) {
        if (*oi == op_next)
            new_operators.push_back(op);
        if (*oi != op)
            new_operators.push_back(*oi);
    }
    // handle move to the end of the operator list
    if (!op_next)
        new_operators.push_back(op);

    m_operators = new_operators;

    emit_plan_changed();
}

void MorphPlan::emit_plan_changed() {
    if (!in_restore) {
        signal_plan_changed();
    }
}

string MorphPlan::generate_id() {
    string chars = id_chars();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, (int)chars.size());

    string id;
    for (size_t i = 0; i < 20; i++)
        id += chars[(size_t)distr(gen)];

    return id;
}

string MorphPlan::id_chars() {
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.:,;/&%$*+@?!#|{}()[]<>=^";
}

Project* MorphPlan::project() {
    return m_project;
}

std::string MorphPlan::id() const {
    return m_id;
}
