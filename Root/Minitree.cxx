#include <FatLeptons/Minitree.h>

using namespace std;

Minitree::Minitree() {
    m_tree = new TTree("tree", "tree");

    makeBranches(m_tree);
}

void Minitree::fillTree(fatLeptonsEvent evt) {
    mj = evt.mj;
    ptj = evt.ptj;
    emfrac = evt.emfrac;
    ntracksj = evt.ntracksj;
    nmatchedels = evt.nmatchedels;
    weight = evt.weight;

    m_tree->Fill();
}

void Minitree::makeBranches(TTree *tree) {
    tree->Branch("mj", &mj, "mj/F");
    tree->Branch("ptj", &ptj, "ptj/F");
    tree->Branch("emfrac", &emfrac, "emfrac/F");
    tree->Branch("ntracksj", &ntracksj, "ntracksj/I");
    tree->Branch("nmatchedels", &nmatchedels, "nmatchedels/I");
    tree->Branch("weight", &weight, "weight/F");
}

void Minitree::setDirectory(TFile *outFile) {
    m_tree->SetDirectory(outFile);
}
