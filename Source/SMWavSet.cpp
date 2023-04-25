// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMWavSet.h"
#include "SMInFile.h"
#include "SMMemOut.h"
#include "SMOutFile.h"
#include "SMUtils.h"

#include <map>
#include <set>

#include <assert.h>

using std::map;
using std::set;
using std::string;
using std::vector;

using namespace SpectMorph;

Error WavSet::save(const string& filename, bool embed_models) {
    OutFile of(filename.c_str(), "SpectMorph::WavSet", SPECTMORPH_BINARY_FILE_VERSION);
    if (!of.open_ok()) {
        fprintf(stderr, "error: can't open output file '%s'.\n", filename.c_str());
        return Error::Code::FILE_NOT_FOUND;
    }

    of.write_string("name", name);
    of.write_string("short_name", short_name);

    for (size_t i = 0; i < waves.size(); i++) {
        of.begin_section("wave");
        of.write_int("midi_note", waves[i].midi_note);
        of.write_int("channel", waves[i].channel);
        of.write_int("velocity_range_min", waves[i].velocity_range_min);
        of.write_int("velocity_range_max", waves[i].velocity_range_max);
        of.write_string("path", waves[i].path);
        if (waves[i].audio) {
            vector<unsigned char> data;

            MemOut mem_out(&data);
            waves[i].audio->save(&mem_out);
            of.write_blob("audio", &data[0], data.size());
        } else if (embed_models) {
            FILE* in = fopen(waves[i].path.c_str(), "rb");
            if (in) {
                vector<unsigned char> data;
                vector<unsigned char> buffer(1024);
                size_t len;

                do {
                    len = fread(&buffer[0], 1, buffer.size(), in);
                    if ((long)len > 0)
                        data.insert(data.end(), buffer.begin(), buffer.begin() + (long)len);
                } while (len > 0);

                if (ferror(in)) {
                    fprintf(stderr, "wavset save: error reading file: %s\n", waves[i].path.c_str());
                } else {
                    of.write_blob("audio", &data[0], data.size());
                }
                fclose(in);
            } else {
                fprintf(stderr, "wavset save: missing file: %s\n", waves[i].path.c_str());
            }
        }
        of.end_section();
    }
    return Error::Code::NONE;
}

Error WavSet::load(const string& filename) {
    clear(); // delete old contents (if any)

    map<string, Audio*> blob_map;

    WavSetWave* wave = nullptr;

    string path;
    string src = "/sources";
    std::string::size_type i = filename.find(src);
    if (i != std::string::npos) {
        string fn = filename;
        path = sm_get_install_dir(INSTALL_DIR_ICLOUD) + fn.erase(i, src.length());
    } else {
        path = sm_get_install_dir(INSTALL_DIR_INSTRUMENTS) + filename;
    }
    InFile ifile(path);
    string section;

    if (!ifile.open_ok())
        return Error::Code::FILE_NOT_FOUND;

    if (ifile.file_type() != "SpectMorph::WavSet")
        return Error::Code::FORMAT_INVALID;

    if (ifile.file_version() != SPECTMORPH_BINARY_FILE_VERSION)
        return Error::Code::FORMAT_INVALID;

    while (ifile.event() != InFile::END_OF_FILE) {
        if (ifile.event() == InFile::BEGIN_SECTION) {
            assert(section == "");
            section = ifile.event_name();

            if (section == "wave") {
                assert(wave == nullptr);
                wave = new WavSetWave();
            }
        } else if (ifile.event() == InFile::END_SECTION) {
            if (section == "wave") {
                assert(wave);

                waves.push_back(*wave);
                delete wave;
                wave = nullptr;
            }

            assert(section != "");
            section = "";
        } else if (ifile.event() == InFile::INT) {
            if (section == "wave") {
                if (ifile.event_name() == "midi_note") {
                    assert(wave);
                    wave->midi_note = ifile.event_int();
                } else if (ifile.event_name() == "channel") {
                    assert(wave);
                    wave->channel = ifile.event_int();
                } else if (ifile.event_name() == "velocity_range_min") {
                    assert(wave);
                    wave->velocity_range_min = ifile.event_int();
                } else if (ifile.event_name() == "velocity_range_max") {
                    assert(wave);
                    wave->velocity_range_max = ifile.event_int();
                } else
                    printf("unhandled int %s %s\n", section.c_str(), ifile.event_name().c_str());
            } else {
                assert(false);
            }
        } else if (ifile.event() == InFile::STRING) {
            if (section == "wave") {
                if (ifile.event_name() == "path") {
                    assert(wave);
                    wave->path = ifile.event_data();
                } else
                    printf("unhandled string %s %s\n", section.c_str(), ifile.event_name().c_str());
            } else if (section == "") {
                if (ifile.event_name() == "name") {
                    name = ifile.event_data();
                } else if (ifile.event_name() == "short_name") {
                    short_name = ifile.event_data();
                } else
                    printf("unhandled global string %s\n", ifile.event_name().c_str());
            } else {
                assert(false);
            }
        } else if (ifile.event() == InFile::BLOB) {
            if (section == "wave") {
                if (ifile.event_name() == "audio") {
                    assert(wave);
                    assert(!wave->audio);

                    GenericIn* blob_in = ifile.open_blob();

                    wave->audio = new Audio();
                    wave->audio->load(blob_in);

                    delete blob_in; // close input file

                    blob_map[ifile.event_blob_sum()] = wave->audio;
                } else {
                    printf("unhandled string %s %s\n", section.c_str(), ifile.event_name().c_str());
                }
            } else {
                assert(false);
            }
        } else if (ifile.event() == InFile::BLOB_REF) {
            if (section == "wave") {
                if (ifile.event_name() == "audio") {
                    assert(wave);
                    assert(!wave->audio);

                    wave->audio = blob_map[ifile.event_blob_sum()];

                    assert(wave->audio);
                } else {
                    printf("unhandled string %s %s\n", section.c_str(), ifile.event_name().c_str());
                }
            } else {
                assert(false);
            }
        } else if (ifile.event() == InFile::READ_ERROR) {
            delete wave;
            return Error::Code::PARSE_ERROR;
        } else {
            assert(false);
        }
        ifile.next_event();
    }
    delete wave;
    return Error::Code::NONE;
}

WavSetWave::WavSetWave() {
    audio = nullptr;
    channel = 0;
    midi_note = -1;
    velocity_range_min = 0;
    velocity_range_max = 127;
}

void WavSet::clear() {
    set<Audio*> to_delete;

    for (vector<WavSetWave>::iterator wi = waves.begin(); wi != waves.end(); wi++) {
        if (wi->audio) {
            to_delete.insert(wi->audio);
            wi->audio = nullptr;
        }
    }

    /* This is a little piece of extra logic that prevents deleting one Audio entry twice;
       the same entry can be reused in one WavSet multiple times (for instance for different
       velocity levels) - and while BLOB_REFs care for reusing the same object in this case,
       we need to ensure here we don't delete the same object twice */
    for (set<Audio*>::const_iterator di = to_delete.begin(); di != to_delete.end(); di++)
        delete *di;

    // now that everything has been delete-d, we can reset the waves vector
    waves.clear();
}

WavSet::~WavSet() {
    clear();
}
