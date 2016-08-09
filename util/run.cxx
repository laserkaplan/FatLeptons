#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/ScanDir.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/LSFDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopAlgs/NTupleSvc.h"
#include "EventLoop/OutputStream.h"

#include "TChain.h"

#include "FatLeptons/FatLeptons.h"

using namespace std;

bool cmdline(int argc, char **argv, map<string,string> &opts);
void usage();

int main(int argc, char **argv) {
    map<string, string> opts;
    if (!cmdline(argc, argv, opts)) return -1;

    if (opts["out"] == "" || opts["in"] == "") {
        cout << "Need in and out!" << endl;
        return -1;
    }

    xAOD::Init().ignore();

    SH::SampleHandler sh;

    std::vector<std::string> ins;
    istringstream instream(opts["in"]);
    string in;
    while (std::getline(instream, in, ',')) {
        ins.push_back(in);
    }

    TChain chain("CollectionTree");
    for (auto &f : ins) chain.Add(f.c_str());
    sh.add(SH::makeFromTChain(opts["out"].c_str(), chain));

    //SH::ScanDir().scan(sh, opts["in"]);

    //SH::scanDQ2(sh, opts["in"]);

    sh.setMetaString("nc_tree", "CollectionTree");

    sh.print();

    EL::Job job;
    job.sampleHandler(sh);

    EL::OutputStream output("minitree");
    job.outputAdd(output);
    EL::NTupleSvc *ntuple = new EL::NTupleSvc("minitree");
    job.algsAdd(ntuple);

    FatLeptons *analysis = new FatLeptons();
    job.algsAdd(analysis);

    analysis->outputName = "minitree";

    EL::DirectDriver driver;
    driver.submit(job, opts["out"]);

    //EL::LSFDriver driver;
    //job.options()->setString(EL::Job::optSubmitFlags, "-q wisc");
    //job.options()->setDouble(EL::Job::optCacheSize, 10*1024*1024);
    //job.options()->setDouble(EL::Job::optCacheLearnEntries, 20);
    //driver.options()->setDouble(EL::Job::optFilesPerWorker, 5);
    //driver.options()->setString(EL::Job::optSubmitFlags, "-L /bin/bash");
    //driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase && source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh";
    //driver.submitOnly(job, opts["out"]);
        
    //EL::PrunDriver driver;
    //driver.options()->setString("nc_outputSampleName", "user.lkaplan." + opts["out"] + ".%in:name[2]");
    //driver.options()->setString("nc_mergeOutput", "true");
    //driver.options()->setDouble("nc_nFilesPerJob", 5);
    //driver.submitOnly(job, opts["out"]);

    return 0;
}

bool cmdline(int argc, char **argv, map<string,string> &opts) {
    opts.clear();

    // defaults
    opts["out"] = "";
    opts["in"] = "";

    for (int i=1;i<argc;i++) {
        string opt=argv[i];

        if (opt=="--help") {usage(); return false;}

        if (0 != opt.find("--")) {
            cout<<"ERROR: options start with '--'!"<<endl;
            return false;
        }
        opt.erase(0,2);
        if (opts.find(opt) == opts.end()) {
            cout<<"ERROR: invalid option '"<<opt<<"'!"<<endl;
            return false;
        }
        string nxtopt=argv[i+1];
        if (0 == nxtopt.find("--") || i+1 >= argc) {
            cout<<"ERROR: option '"<<opt<<"' requires value!"<<endl;
            return false;
        }

        opts[opt] = nxtopt;
        i++;
    }

    return true;
}

void usage() {
    cout<<"USAGE: run [-option value]\n\n"
        <<"options [default]:\n\n"
        <<"--help\n"
        <<"--out (required!)\n"
        <<"--in (required!)\n"
        <<endl;
}
