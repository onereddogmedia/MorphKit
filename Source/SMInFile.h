// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_INFILE_HH
#define SPECTMORPH_INFILE_HH

#include <set>
#include <stdio.h>
#include <string>
#include <vector>

#include "SMAudio.h"
#include "SMMMapIn.h"
#include "SMStdioIn.h"

namespace SpectMorph {

/**
 * \brief Class to read SpectMorph binary data.
 *
 * This class allows reading SpectMorph binary data from a file (with name) or
 * from a GenericIn object; the files consist of events, which are read one by
 * one until an error occurs or until EOF. The events are typed, end the event
 * data should be queried according to the current event type.
 */
class InFile {
  public:
    enum Event {
        NONE,
        END_OF_FILE,
        READ_ERROR,
        BEGIN_SECTION,
        END_SECTION,
        BOOL,
        INT,
        STRING,
        FLOAT,
        FLOAT_BLOCK, // not used
        UINT16_BLOCK,
        BLOB,
        BLOB_REF
    };

  protected:
    GenericIn* file;
    bool file_delete;
    Event current_event;
    std::string current_event_str;
    bool current_event_bool;
    int current_event_int;
    std::string current_event_data;
    float current_event_float;
    std::vector<float> current_event_float_block;
    std::vector<uint16_t> current_event_uint16_block;
    size_t current_event_blob_pos;
    size_t current_event_blob_size;
    std::string current_event_blob_sum;
    std::string m_file_type;
    int m_file_version;

    bool read_raw_bool(bool& b);
    bool read_raw_string(std::string& str);
    bool read_raw_int(int& i);
    bool read_raw_float(float& f);
    bool read_raw_float_block(std::vector<float>& fb);
    bool read_raw_uint16_block(std::vector<uint16_t>& ib);

    void read_file_type_and_version();

  public:
    InFile(const std::string& filename);
    InFile(GenericIn* file);
    ~InFile();

    /**
     * Check if file open succeeded.
     *
     * \returns true if file was opened successfully, false otherwise
     */
    bool open_ok() {
        return file != nullptr;
    }
    Event event();
    std::string event_name();
    float event_float();
    int event_int();
    bool event_bool();
    std::string event_data();
    const std::vector<float>& event_float_block();
    const std::vector<uint16_t>& event_uint16_block();
    std::string event_blob_sum();

    void next_event();
    std::string file_type();
    int file_version();

    GenericIn* open_blob();
};

} // namespace SpectMorph

#endif /* SPECTMORPH_INFILE_HH */
