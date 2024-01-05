// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_OUTPUT_HH
#define SPECTMORPH_MORPH_OUTPUT_HH

#include "SMMath.h"
#include "SMModulationList.h"
#include "SMMorphOperator.h"
#include "SMProperty.h"
#include "SMUtils.h"

#include <string>

namespace SpectMorph {

class MorphOutput;

class MorphOutput : public MorphOperator {
  public:
    struct Config : public MorphOperatorConfig {
        std::vector<MorphOperatorPtr> channel_ops;
    };
    Config m_config;

  protected:
    std::vector<std::string> load_channel_op_names;

  public:
    MorphOutput(MorphPlan* morph_plan);
    ~MorphOutput() override;

    // inherited from MorphOperator
    const char* type() override;
    int insert_order() override;
    OutputType output_type() override;

    std::vector<MorphOperator*> dependencies() override;
    MorphOperatorConfig* clone_config() override;

    void set_channel_op(int ch, MorphOperator* op);
    MorphOperator* channel_op(int ch);

    /* slots: */
    void on_operator_removed(MorphOperator* op);
};

} // namespace SpectMorph

#endif
