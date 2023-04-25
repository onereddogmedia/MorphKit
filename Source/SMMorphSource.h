// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_SOURCE_HH
#define SPECTMORPH_MORPH_SOURCE_HH

#include "SMMorphOperator.h"

#include <string>

namespace SpectMorph {

class MorphSource : public MorphOperator {
  public:
    struct Config : public MorphOperatorConfig {
        std::string path;
    };
    Config m_config;

  protected:
    std::string m_smset;

  public:
    MorphSource(MorphPlan* morph_plan);
    ~MorphSource() override;

    // inherited from MorphOperator
    const char* type() override;
    int insert_order() override;
    OutputType output_type() override;
    MorphOperatorConfig* clone_config() override;

    void set_smset(const std::string& smset);
    std::string smset();
};

} // namespace SpectMorph

#endif
