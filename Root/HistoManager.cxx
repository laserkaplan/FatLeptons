#include "FatLeptons/HistoManager.h"

using namespace std;

HistoManager::HistoManager() {
}

HistoManager::~HistoManager() {
    for (auto a : histos) {
        delete a;
    }
}
