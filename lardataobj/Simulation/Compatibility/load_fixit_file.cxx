#include "TFile.h"

#include "cetlib/search_path.h"

namespace {
  int load_fixit_file()
  {
    cet::search_path const sp{"FW_SEARCH_PATH", std::nothrow};
    if (std::string filename; sp.find_file("simIDE_streamer_info.root", filename)) {
      TDirectory::TContext cx;
      TFile file{filename.c_str()};
      return 0;
    }
    return 1;
  }

  auto rc = load_fixit_file();
}
