// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMProperty.h"
#include "SMModulationList.h"
#include "SMMorphPlan.h"

using namespace SpectMorph;

using std::string;

Property::Property(MorphOperator* op, const string& identifier) : m_op(op), m_identifier(identifier) {
    m_op->register_property(this);
}

Property::~Property() {
}

void Property::set_modulation_data(ModulationData* mod_data) {
    Range r = float_range();
    mod_data->min_value = (float)r.min_value;
    mod_data->max_value = (float)r.max_value;

    auto property_scale = float_scale();
    mod_data->property_scale = property_scale;

    switch (property_scale) {
        case Scale::LOG:
            mod_data->value_scale = (float)log2f((float)(r.max_value / r.min_value));
            break;
        case Scale::LINEAR:
            mod_data->value_scale = (float)(r.max_value - r.min_value);
            break;
        case Scale::NONE:
        default:
            mod_data->value_scale = 0;
            break;
    }

    m_modulation_list = std::make_unique<ModulationList>(*mod_data, *this);

    connect(m_modulation_list->signal_modulation_changed, [this]() { signal_modulation_changed(); });
}
