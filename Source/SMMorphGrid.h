// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_GRID_HH
#define SPECTMORPH_MORPH_GRID_HH

#include "SMMorphOperator.h"

#include <map>

namespace SpectMorph {

struct MorphGridNode {
    std::string smset; // a node has an instrument (smset)
    std::string path;

    MorphGridNode();
};

class MorphGrid : public MorphOperator {
  public:
    struct Config : public MorphOperatorConfig {
        unsigned int width;
        unsigned int height;

        double x_morphing;
        double y_morphing;

        ControlType x_control_type;
        ControlType y_control_type;
        ControlType node_a_db_control_type;
        ControlType node_b_db_control_type;
        ControlType node_c_db_control_type;
        ControlType node_d_db_control_type;

        std::vector<std::vector<MorphGridNode>> input_node;
    };

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

    void set_width(unsigned int width);
    unsigned int width();

    void set_height(unsigned int height);
    unsigned int height();

    void set_x_morphing(double new_value);
    void set_y_morphing(double new_value);
    void set_x_control_type(ControlType new_control_type);
    void set_y_control_type(ControlType new_control_type);
    void set_node_a_db_control_type(ControlType new_control_type);
    void set_node_b_db_control_type(ControlType new_control_type);
    void set_node_c_db_control_type(ControlType new_control_type);
    void set_node_d_db_control_type(ControlType new_control_type);

    void set_input_node(unsigned int x, unsigned int y, const MorphGridNode& node);
    MorphGridNode input_node(unsigned int x, unsigned int y);
};

} // namespace SpectMorph

#endif
