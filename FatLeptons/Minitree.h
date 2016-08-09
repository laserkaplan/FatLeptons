#ifndef FatLeptons_Minitree_H
#define FatLeptons_Minitree_H

#include "TFile.h"
#include "TTree.h"

struct fatLeptonsEvent {
    float mj;
    float ptj;
    float emfrac;
    int ntracksj;
    int nmatchedels;
    float weight;
};

class Minitree {
    public:
        Minitree();
        
        void fillTree(fatLeptonsEvent evt);
        void setDirectory(TFile *outFile);

    private:
        void makeBranches(TTree *tree);

        // tree
        TTree *m_tree; //!

        // branches
        float mj; //!
        float ptj; //!
        float emfrac; //!
        int ntracksj; //!
        int nmatchedels; //!
        float weight; //!
};

#endif
