#pragma once
#include <iostream>
typedef int64_t LL;

class Progress_printer{

    public:

    LL n_jobs;
    LL processed;
    LL total_prints;
    LL next_print;

    Progress_printer(LL n_jobs, LL total_prints) : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0) {}

    void job_done(string msg = ""){
        if(next_print == processed){
            LL progress_percent = round(100 * ((double)processed / n_jobs));
            std::cerr << to_string(progress_percent) + "%" << " " << msg << std::endl;
            next_print += n_jobs / total_prints;
        }
        processed++;
    }

};