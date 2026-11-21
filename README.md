# Dimensijos mažinimas klasterizavime

## Tikslas: Išskirti ir apibrėžti deimantų duomenų aibėje susidariusius klasterius

### Darbo plano aprašymas:

1. Paruošti duomenų aibę
    - Užpildyti praleistas reikšmes
    - Panaikinti išskirtis
    - Panaikinti koreliuotus požymius
    - Sunormuoti
1. Sukurti 5 aibes tolesniam darbui:
    - Orginali 
    - PCA (5-6 dimencijų)
    - PCA (2 dimencijų)
    - UMAP (5-6 dimencijų)
    - UMAP (2 dimencijų)
1. Pasirinkti požymius, pagal kuriuos bus atliekamas klasterizavimas (atrinkti požymiai kurių kor. kof. < 95%)
1. Nustatyti optimalų klasterių skaičių naudojant:
    - Empyrinį metodą
    - Elbow metodą
    - Vidutinio silueto
    - Kelio metodą
    - Dendogramą
1. Kiekvienu klasifikavimo algoritmu (K-means, hierarchinis, DBSCAN) praeiti pro visas aibes ir:
    - Apskaičiuoti:
        - Vidinius matus:
            -  Silhouette
            -  Calinski–Harabasz
            -  Davies–Bouldin
        - Išorinius matus:
            - Jaccard
            - ARI
            - NMI
        - Stabilumo eksperimentas    
    - Bendrai visiems:
        - Sudaryti graifkų matricą
        - Įrašyti matus į lentelę
1. Interpretuoti ir aprašyti patį geriausią algoritmą 
