// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_OPERATOR_HH
#define SPECTMORPH_MORPH_OPERATOR_HH

#include "SMInFile.h"
#include "SMOutFile.h"
#include "SMSignal.h"
#include "SMWavSetRepo.h"

#include <map>
#include <memory>

namespace SpectMorph {

struct MorphOperatorConfig {
    MorphOperatorConfig() : wav_set_repo(nullptr) {
    }
    MorphOperatorConfig(const MorphOperatorConfig& other) : wav_set_repo(other.wav_set_repo) {
    }
    MorphOperatorConfig& operator=(MorphOperatorConfig&& other) {
        if (this != &other) {
            wav_set_repo = other.wav_set_repo;
        }
        return *this;
    }
    virtual ~MorphOperatorConfig();
    WavSetRepo* wav_set_repo;
};

typedef std::shared_ptr<MorphOperatorConfig> MorphOperatorConfigP;

class MorphOperatorView;
class MorphPlan;
class MorphOperatorPtr;

class MorphOperator : public SignalReceiver {
  public:
    typedef std::map<std::string, MorphOperator*> OpNameMap;

  protected:
    MorphPlan* m_morph_plan;
    std::string m_name;
    std::string m_id;

  public:
    enum OutputType { OUTPUT_NONE, OUTPUT_AUDIO, OUTPUT_CONTROL };
    enum ControlType {
        CONTROL_GUI = 1,
        CONTROL_SIGNAL_1 = 2,
        CONTROL_SIGNAL_2 = 3,
        CONTROL_OP = 4,

        /* note: don't reorder items here, as we need to be compatible with old files */
        CONTROL_SIGNAL_3 = 5,
        CONTROL_SIGNAL_4 = 6,
        CONTROL_SIGNAL_5 = 7,
        CONTROL_SIGNAL_6 = 8
    };
    typedef uintptr_t PtrID;

    MorphOperator(MorphPlan* morph_plan);
    virtual ~MorphOperator();

    virtual const char* type() = 0;
    virtual int insert_order() = 0;
    virtual OutputType output_type() = 0;
    virtual std::vector<MorphOperator*> dependencies();
    virtual MorphOperatorConfig* clone_config() = 0;

    MorphPlan* morph_plan();

    std::string type_name();

    std::string name();
    void set_name(const std::string& name);

    bool can_rename(const std::string& name);

    std::string id();
    void set_id(const std::string& id);

    PtrID ptr_id() const {
        /* ptr_id is derived from MorphOperator*, which means that for a given
         * MorphPlan, these ids never collide, but if the plan is modified,
         * the same ptr_id can be taken by a new MorphOperator*
         */
        return PtrID(this);
    }

    static MorphOperator* create(const std::string& type, MorphPlan* plan);
};

class MorphOperatorPtr {
  private:
    MorphOperator* m_ptr = nullptr;

  public:
    operator bool() const {
        return m_ptr != nullptr;
    }
    MorphOperator* get() const {
        return m_ptr;
    }

    MorphOperator::PtrID ptr_id() const {
        return MorphOperator::PtrID(m_ptr);
    }

    void set(MorphOperator* ptr) {
        m_ptr = ptr;
    }
};

} // namespace SpectMorph

#endif
