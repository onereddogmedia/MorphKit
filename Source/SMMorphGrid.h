// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_GRID_HH
#define SPECTMORPH_MORPH_GRID_HH

#include "SMMorphOperator.h"

#include <map>

namespace SpectMorph {

struct MorphGridNode {
    MorphOperatorPtr op; // a node has either an operator (op) as input,
    std::string smset;   // or an instrument (smset)
    WavSet* wav_set = nullptr;

    MorphGridNode();
};

class MorphGrid : public MorphOperator {
  public:
    struct Config : public MorphOperatorConfig {
        int width;
        int height;

        double x_morphing;
        double y_morphing;

        std::vector<std::vector<MorphGridNode>> input_node;
    };
    static constexpr auto P_X_MORPHING = "x_morphing";
    static constexpr auto P_Y_MORPHING = "y_morphing";

  protected:
    Config m_config;

    std::map<std::string, std::string> load_map;

    void update_size();

  public:
    MorphGrid(MorphPlan* morph_plan);
    ~MorphGrid() override;

    // inherited from MorphOperator
    const char* type() override;
    int insert_order() override;
    OutputType output_type() override;
    MorphOperatorConfig* clone_config() override;

    std::vector<MorphOperator*> dependencies() override;

    void set_width(int width);
    int width();

    void set_height(int height);
    int height();

    double x_morphing();
    void set_x_morphing(double new_value);
    double y_morphing();
    void set_y_morphing(double new_value);

    void set_input_node(int x, int y, const MorphGridNode& node);
    MorphGridNode input_node(int x, int y);
    std::string input_node_label(int x, int y);

    /* slots: */
    void on_operator_removed(MorphOperator* op);
};

} // namespace SpectMorph

#endif
