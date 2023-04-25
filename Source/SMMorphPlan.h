// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_PLAN_HH
#define SPECTMORPH_MORPH_PLAN_HH

#include "SMAudio.h"
#include "SMMorphOperator.h"
#include "SMObject.h"
#include "SMSignal.h"
#include "SMUtils.h"

namespace SpectMorph {

class Project;

class MorphPlan : public Object {
  public:
    class ExtraParameters {
      public:
        virtual ~ExtraParameters() = default;
        virtual std::string section() = 0;
        virtual void save(OutFile& out_file) = 0;
        virtual void handle_event(InFile& ifile) = 0;
    };

    constexpr static int N_CONTROL_INPUTS = 8;

  protected:
    Project* m_project = nullptr;
    std::vector<MorphOperator*> m_operators;
    std::string m_id;

    bool in_restore;

  public:
    MorphPlan(Project& project);
    ~MorphPlan() override;

    Project* project();
    std::string id();

    enum AddPos { ADD_POS_AUTO, ADD_POS_END };

    void add_operator(MorphOperator* op, AddPos = ADD_POS_END, const std::string& name = "",
                      const std::string& id = "");
    const std::vector<MorphOperator*>& operators();
    void clear();
    void remove(MorphOperator* op);
    void move(MorphOperator* op, MorphOperator* op_next);

    void emit_plan_changed();

    static std::string id_chars();
    static std::string generate_id();

    Signal<> signal_plan_changed;
    Signal<> signal_need_view_rebuild;
    Signal<MorphOperator*> signal_operator_removed;
    Signal<MorphOperator*> signal_operator_added;
};

typedef RefPtr<MorphPlan> MorphPlanPtr;

} // namespace SpectMorph

#endif
