// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_PLAN_HH
#define SPECTMORPH_MORPH_PLAN_HH

#include "SMAudio.h"
#include "SMMorphOperator.h"
#include "SMSignal.h"
#include "SMUtils.h"

namespace SpectMorph {

class Project;

class MorphPlan : public SignalReceiver {
  public:
    class ExtraParameters {
      public:
        virtual ~ExtraParameters() = default;
        virtual std::string section() = 0;
        virtual void handle_event(InFile& ifile) = 0;
    };

    constexpr static int N_CONTROL_INPUTS = 8;

  protected:
    Project* m_project = nullptr;
    std::vector<MorphOperator*> m_operators;
    std::string m_id;

    WavSetRepo* m_wav_set_repo;

    bool in_restore;

    void clear();
    Error load_internal(GenericIn* in, ExtraParameters* params = nullptr);

  public:
    MorphPlan(Project& project);
    ~MorphPlan() override;

    Project* project();
    std::string id() const;

    enum AddPos { ADD_POS_AUTO, ADD_POS_END };

    void add_operator(MorphOperator* op, AddPos = ADD_POS_END, const std::string& name = "", const std::string& id = "",
                      bool load_folded = false);
    const std::vector<MorphOperator*>& operators() const;
    void remove(MorphOperator* op);
    void move(MorphOperator* op, MorphOperator* op_next);

    void emit_plan_changed();

    void reloadWavSet();

    static std::string id_chars();
    static std::string generate_id();

    WavSetRepo* wavSetRepo() {
        return m_wav_set_repo;
    }

    Signal<> signal_plan_changed;
    Signal<> signal_index_changed;
    Signal<> signal_need_view_rebuild;
    Signal<MorphOperator*> signal_operator_removed;
    Signal<MorphOperator*> signal_operator_added;
};

} // namespace SpectMorph

#endif
