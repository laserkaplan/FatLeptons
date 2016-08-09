#ifndef HistoManager_H
#define HistoManager_H

#include <map>
#include <vector>

#include "TH1.h"

class HistoManager {
    public:
        HistoManager();
        virtual ~HistoManager();

        template <class T>
        void book(T hist, std::string name) {
            histos.push_back(hist);
            nameMap.insert(std::pair<std::string, int>(name, (int) histos.size() - 1));
        }
        template <typename T>
        void fill(std::string name, T val) {
            histos[nameMap[name]]->Fill(val, m_weight);
        }
        std::vector<TH1*> getHistos() {return histos;};
        void setWeight(float weight) {m_weight = weight;};

    protected:
        std::vector<TH1*> histos;
        std::map<std::string, int> nameMap;
        float m_weight;
};

#endif
