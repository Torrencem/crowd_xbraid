
#include "tribraid.hpp"
#include <type_traits>
#include "split_line_search.c"

class TriBraidLineSearchApp : public TriBraidApp {
public:
    TriBraidLineSearchApp(MPI_Comm comm_t_, double tstart_, double tstop_, int ntime_)
        : TriBraidApp(comm_t_, tstart_, tstop_, ntime_) {}

    virtual int Sync(BraidSyncStatus &sstatus_) override {
        const bool is_standard_layout_v = std::is_standard_layout<BraidSyncStatus>::value;
        static_assert(is_standard_layout_v, "BraidSyncStatus is not standard layout, the following cast is UB:");
        braid_SyncStatus *sstatus = reinterpret_cast<braid_SyncStatus *>(&sstatus_);
        braid_App app = (braid_App) this;
        return line_search_sync(app, *sstatus);
    }
};
