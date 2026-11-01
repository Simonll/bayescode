// AcceptanceTracker: generic MCMC acceptance-rate reporter.
//
// Usage: enable with env var TRACE_ACCEPTANCE=1. Each move calls
//   AcceptanceTracker::get().add(label, nacc, ntot);
// and, once per MCMC iteration, calls
//   AcceptanceTracker::get().periodicReport(iteration);
// which prints cumulative acceptance per label every `period` iterations.
//
// No dependency on spdlog/logging.hpp (kept header-only and self-contained),
// consistent with the models' existing use of std::cout for diagnostics.
#pragma once

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

class AcceptanceTracker {
  public:
    static AcceptanceTracker& get() {
        static AcceptanceTracker instance;
        return instance;
    }

    // true if env TRACE_ACCEPTANCE is set and non-empty
    bool enabled() const { return std::getenv("TRACE_ACCEPTANCE") != nullptr; }

    void setPeriod(int p) { period_ = p; }

    // accumulate a move's outcomes under a label; no-op if disabled
    void add(const std::string& label, double nacc, double ntot) {
        if (!enabled()) { return; }
        nacc_[label] += nacc;
        ntot_[label] += ntot;
    }

    // print cumulative acceptance for all labels (and reset) if this call is
    // a multiple of period_ from the last report. Accepts the current external
    // iteration counter for the header line.
    void periodicReport(int iter, bool force = false) {
        if (!enabled()) { return; }
        if (period_ <= 0) { period_ = 100; }
        if (!force && (iter - last_iter_ < period_)) { return; }
        last_iter_ = iter;
        for (const auto& kv : ntot_) {
            const std::string& label = kv.first;
            double ntot = kv.second;
            double nacc = nacc_.count(label) ? nacc_[label] : 0.0;
            double rate = ntot > 0 ? nacc / ntot : 0.0;
            std::cout << "[accept] " << label << ": " << rate
                      << " (" << static_cast<long>(nacc) << "/"
                      << static_cast<long>(ntot) << ")\n";
        }
        nacc_.clear();
        ntot_.clear();
    }

  private:
    AcceptanceTracker() {
        // default period (overridable via setPeriod)
    }
    ~AcceptanceTracker() = default;
    AcceptanceTracker(const AcceptanceTracker&) = delete;
    AcceptanceTracker& operator=(const AcceptanceTracker&) = delete;

    int period_{100};
    int last_iter_{0};
    std::map<std::string, double> nacc_;
    std::map<std::string, double> ntot_;
};
