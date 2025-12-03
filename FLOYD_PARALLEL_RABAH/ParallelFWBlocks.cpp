#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <iostream>
#include "ParallelFWBlocks.hpp"
#include "Distribution.hpp"

using namespace std;

static const int INF = 1000000000;

// ===== HELPERS  =====

// Cette fonction c'est moi qui l'utilise pour faire le Floyd-Warshall
// mais seulement a l'intérieur d'un seul bloc (un carre b x b).
// bs c'est la "vraie" taille du bloc, parce que des fois a la fin
// le bloc déborde de la matrice et du coup il est pas complet.

// En gros ici je fais les 3 boucles classiques de Floyd-Warshall :
// kk c'est le pivot dans le bloc,
// i c'est la ligne,
// j c'est la colonne.

// Pour chaque i et j, je regarde si passer par kk (donc faire i -> kk -> j)
// donne un chemin plus court que ce que j'avais deja.

// Je fais bien gaffe aux INF : si la distance i->kk ou kk->j vaut INF,
// ça veut dire que ce chemin existe pas, donc je continue et j’évite
// de faire des additions qui servent a rien.

// Au final, si "via" est plus petit que la distance actuelle i->j,
// alors je mets à jour la case dans Dkk.
static void fw_block(int* Dkk, int bs, int b) {
    for (int kk = 0; kk < bs; ++kk) {
        for (int i = 0; i < bs; ++i) {
            int dik = Dkk[i * b + kk];
            if (dik == INF) continue;
            for (int j = 0; j < bs; ++j) {
                int kkj = Dkk[kk * b + j];
                if (kkj == INF) continue;
                int via = dik + kkj;
                int& dij = Dkk[i * b + j];
                if (via < dij) dij = via;
            }
        }
    }
}

// Ici je fais la mise à jour des blocs qui sont sur la même ligne
// que le bloc pivot Dkk. Donc en gros : [ D(k,k) ] -> propage vers [ D(k,J) ].
//
// Cette fonction fw_row prend le bloc pivot Dkk et un bloc DkJ
// qui est dans la même "row" (ligne) que le pivot.
//
// bs c'est la taille réelle du bloc pivot (peut être plus petit en bas à droite),
// et wJ c'est la largeur réelle du bloc DkJ (pareil, pour éviter de dépasser la vraie matrice).
//
// Dans cette boucle, je fais exactement le même principe que FW classique,
// sauf que cette fois j'utilise les distances du pivot Dkk et j'essaie
// d'améliorer les distances DkJ (donc entre les sommets de la ligne k et la colonne J).
//
// Pour chaque i dans le bloc pivot,
// je teste tous les kk comme pivot local,
// puis tous les j de la bande DkJ.
//
// Si le chemin i -> kk est INF ou kk -> j est INF,
// je skip parce que ça sert à rien de continuer.
// Sinon je calcule le chemin "via" = i -> kk -> j,
// et si il est plus petit que ce que j'avais dans DkJ,
// je mets à jour.
//
// En gros, cette fonction sert juste à propager la ligne du pivot
// vers les blocs de droite dans la grille de blocs.
static void fw_row(const int* Dkk, int* DkJ, int bs, int wJ, int b) {
    for (int i = 0; i < bs; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int dik = Dkk[i * b + kk];
            if (dik == INF) continue;
            for (int j = 0; j < wJ; ++j) {
                int kkj = DkJ[kk * b + j];
                if (kkj == INF) continue;
                int via = dik + kkj;
                int& dij = DkJ[i * b + j];
                if (via < dij) dij = via;
            }
        }
    }
}


// Cette fonction fw_col c'est un peu la même idée que fw_row,
// sauf que là je travaille sur la COLONNE du pivot au lieu de la ligne.
//
// Le bloc Dkk c'est toujours le bloc pivot (au centre),
// mais cette fois je propage ses infos vers un bloc qui est juste
// en dessous du pivot, qu'on appelle Dik.
//
// hI c'est la hauteur réelle du bloc Dik (ça évite de dépasser la matrice
// quand on tombe sur les blocs du bas), et bs c'est la taille réelle du pivot.
//
// Le principe est simple : pour chaque cellule (i,j) du bloc Dik,
// j'essaie de l'améliorer en passant par le pivot : i -> kk -> j.
//
// i = ligne du bloc Dik,
// kk = pivot local dans le bloc Dkk,
// j = colonne dans le bloc Dik.
//
// Comme d'habitude, si ik ou kj vaut INF,
// je continue parce que ça veut dire que ce chemin n'existe pas,
// donc aucun intérêt de calculer.
//
// Sinon je calcule "via" = ik + kj.
// Si "via" est plus petit que Dik[i][j],
// je mets à jour.
//
// Au final fw_col met à jour tous les blocs situés SOUS le bloc pivot,
// en utilisant les infos du pivot pour améliorer leurs distances.
static void fw_col(int* Dik, const int* Dkk, int hI, int bs, int b) {
    for (int i = 0; i < hI; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int ik = Dik[i * b + kk];
            if (ik == INF) continue;
            for (int j = 0; j < bs; ++j) {
                int kj = Dkk[kk * b + j];
                if (kj == INF) continue;
                int via = ik + kj;
                int& ij = Dik[i * b + j];
                if (via < ij) ij = via;
            }
        }
    }
}


// fw_inner c'est la partie où je mets à jour les blocs "au milieu",
// c’est-à-dire ceux qui ne sont ni dans la même ligne que le pivot,
// ni dans la même colonne, mais qui sont dans le rectangle autour.
//
// En fait pour un bloc Dij, je peux l’améliorer en combinant :
//    Dik  (bloc dans la colonne du pivot)
//    DkJ  (bloc dans la ligne du pivot)
//
// Donc je fais le chemin i -> kk via le bloc Dik,
// puis kk -> j via DkJ.
//
// hI c’est la hauteur réelle du bloc (si on est en bas de la matrice),
// wJ c’est la largeur réelle du bloc (si on est à droite),
// et bs c’est la taille du pivot local.
//
// Comme d’hab, si ik est INF ou kj est INF,
// je continue direct, parce que ça veut dire que ce chemin n’existe même pas.
//
// Je calcule "via" = ik + kj,
// et si ça améliore la distance qu’il y avait déjà dans Dij,
// je mets à jour.
//
// En vrai fw_inner c’est la partie qui propage le pivot dans tous les blocs
// qui ne sont pas directement collés au pivot, un peu comme une mise à jour
// du carré central dans Floyd-Warshall mais en version découpée en blocs.
static void fw_inner(const int* Dik, const int* DkJ, int* Dij, 
                     int hI, int wJ, int bs, int b) {
    for (int i = 0; i < hI; ++i) {
        for (int kk = 0; kk < bs; ++kk) {
            int ik = Dik[i * b + kk];
            if (ik == INF) continue;
            for (int j = 0; j < wJ; ++j) {
                int kj = DkJ[kk * b + j];
                if (kj == INF) continue;
                int via = ik + kj;
                int& ij = Dij[i * b + j];
                if (via < ij) ij = via;
            }
        }
    }
}

// ===== =====
int* ParallelFloydWarshallBlocks(int n, int* mat) {
    using namespace std;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // ===== CHOIX DE b : PRIORITÉ À b = n / sqrt(p) SI POSSIBLE =====

      // Ici je choisis la taille des blocs b.
    // L'idée de base, c'est que si j'ai p processus et que p est un carré parfait,
    // je préfère avoir une grille sqrt(p) x sqrt(p) avec des blocs bien réguliers.
    //
    // Du coup je commence par calculer sqrt(p) et je regarde si :
    //  - sqrtp * sqrtp == size (donc p est un carré parfait),
    //  - et en plus n est multiple de sqrtp (pour aligner pile les blocs).
    int sqrtp = (int)std::lround(std::sqrt((double)size));
    bool grilleCarree = (sqrtp * sqrtp == size) && (sqrtp > 0) && (n % sqrtp == 0);

    int b;
    if (grilleCarree) {
         // Cas “propre” : je fais comme dans le cours,
        // chaque dimension de la matrice est découpée en sqrt(p) blocs.
        b = n / sqrtp;
    } else {
        // Sinon, je prends un b “adaptatif”.
        // Je pars de n / sqrt(p) mais en mode ceil, donc j'arrondis vers le haut.
        int denom = std::max(1, sqrtp);
        b = (n + denom - 1) / denom;            // ceil division entière
        b = std::max(32, std::min(256, b));     // évite des blocs trop petits ou trop grands
    }

        // nb = nombre de blocs par dimension (en arrondissant vers le haut si ça ne tombe pas juste)
    int nb = (n + b - 1) / b;

    if (rank == 0) {
        cout << "[INFO] Taille matrice : " << n << "x" << n << endl;
        cout << "[INFO] Taille bloc    : " << b << "x" << b << endl;
        cout << "[INFO] Nombre blocs   : " << nb << "x" << nb << endl;
        cout << "[INFO] Processus      : " << size << endl;
        if (!grilleCarree) {
            cout << "[WARN] p non carré parfait ou n non multiple de sqrt(p) : b adaptatif utilise." << endl;
        }
    }

    // ===== Grille de processus =====

     // Ici je construis une grille 2D de processus (Pr x Pc).
    // Si j'ai la chance d'avoir p carré parfait et n multiple de sqrt(p),
    // je force un hint vers une grille carrée sqrtp x sqrtp.
    // Sinon je laisse MPI_Dims_create choisir automatiquement une décomposition.

    int dims[2] = {0, 0};
    if (grilleCarree) {
        dims[0] = sqrtp;
        dims[1] = sqrtp;
    }
    MPI_Dims_create(size, 2, dims);
    int Pr = dims[0];
    int Pc = dims[1];

    // ===== Distribution des blocs =====

    // Ici je demande quels blocs appartiennent à CE processus.
    // computeLocalBlocks va parcourir tous les blocs (bi,bj)
    // et garder ceux dont le ownerOf(bi,bj,Pr,Pc) == rank.
    vector<BlockInfo> localBlocks = computeLocalBlocks(n, b, Pr, Pc, rank);
    int numLocal = (int)localBlocks.size();
    int blockArea = b * b;

    // Je stocke les données de tous mes blocs locaux dans un gros tableau 1D.
    // Chaque bloc fait b*b cases, donc au total numLocal * blockArea.
    vector<int> localData(numLocal * blockArea);

    // Ce tableau me donne, pour chaque bloc global (bi,bj),
    // l'indice du bloc local correspondant chez CE processus.
    // Si une case vaut -1, ça veut dire que ce processus ne possède pas ce bloc.
    vector<int> localIndex(nb * nb, -1);
    
    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        localIndex[info.bi * nb + info.bj] = idx;
    }

    // ===== Initialisation =====

    // Maintenant je remplis chaque bloc local à partir de la matrice globale mat.
    // Je fais attention aux blocs du bord :
    //  - si gi ou gj dépasse n, je mets INF (padding),
    //  - si on est sur la diagonale (gi == gj), je mets 0,
    //  - si mat[gi*n + gj] == 0, ça veut dire pas d'arête -> INF,
    //  - sinon je recopie le poids de l'arête.


    for (int idx = 0; idx < numLocal; ++idx) {
        const BlockInfo& info = localBlocks[idx];
        int i0 = info.offset_i;
        int j0 = info.offset_j;
        int* blk = &localData[idx * blockArea];

        for (int ii = 0; ii < b; ++ii) {
            int gi = i0 + ii;    // indice global en ligne 
            for (int jj = 0; jj < b; ++jj) {
                int gj = j0 + jj;   // indice global en colonne 
                int val;
                if (gi >= n || gj >= n) {
// On est en dehors de la vraie matrice -> padding INF
                    val = INF;
                } else if (gi == gj) {
 // Distance de i à i = 0

                    val = 0;
                } else if (mat[gi * n + gj] == 0) {
// Pas d'arête explicite -> pas de chemin direct

                    val = INF;
                } else {
 // Arête présente dans la matrice d'adjacence -> on copie le poids

                    val = mat[gi * n + gj];
                }
                blk[ii * b + jj] = val;
            }
        }
    }

    // =====COMMUNICATIONS NON-BLOQUANTES =====
    // Ici je prépare deux gros tableaux pour stocker, pour chaque bloc de la ligne k
// et de la colonne k, les données dont j'ai besoin pour mettre à jour le reste.
//
// rowBlocks[jb] : contient le bloc D(k, jb)
// colBlocks[ib] : contient le bloc D(ib, k)
    vector<vector<int>> rowBlocks(nb, vector<int>(blockArea));
    vector<vector<int>> colBlocks(nb, vector<int>(blockArea));
    
    // Buffers pour les communications asynchrones
    vector<MPI_Request> requests;
    requests.reserve(2 * nb);

  // ===== Boucle principale sur les blocs pivots =====
// kk parcourt les blocs diagonaux (k,k) en coordonnées de blocs.
    for (int kk = 0; kk < nb; ++kk) {
            // Je récupère le rang MPI qui possède le bloc pivot (kk,kk)
        int pivotOwner = ownerOf(kk, kk, Pr, Pc);
        int pivotLocalIdx = localIndex[kk * nb + kk];

            // Je prépare un buffer pour le bloc pivot, initialisé à INF par défaut
        vector<int> pivotBlock(blockArea, INF);
        int bs = std::min(b, n - kk * b);

        // ===== Phase A : Bloc pivot =====
            // Seul le processus qui possède le bloc (kk,kk) fait le Floyd-Warshall local dessus.

        if (rank == pivotOwner && pivotLocalIdx != -1) {
            int* DkkLocal = &localData[pivotLocalIdx * blockArea];
                    // Je fais FW sur le bloc pivot uniquement (fw_block)

            fw_block(DkkLocal, bs, b);
                    // Puis je copie le bloc pivot mis à jour dans pivotBlock

            std::copy(DkkLocal, DkkLocal + blockArea, pivotBlock.begin());
        }
            // Ensuite je diffuse le bloc pivot à tout le monde (broadcast classique).

        MPI_Bcast(pivotBlock.data(), blockArea, MPI_INT, pivotOwner, MPI_COMM_WORLD);
            // Et je le garde aussi comme "k-ième" bloc de ligne et de colonne

        rowBlocks[kk] = pivotBlock;
        colBlocks[kk] = pivotBlock;

        // Phase B.1 : mise à jour de la LIGNE de blocs (k, jb)
    // Idée : le processus qui possède D(k,jb) le met à jour localement,
    // puis on Ibcast ce bloc à tous les autres.
        requests.clear();
        
        for (int jb = 0; jb < nb; ++jb) {
            if (jb == kk) continue;  // on saute le pivot lui-même

            int ownerRow = ownerOf(kk, jb, Pr, Pc);
            int wJ = std::min(b, n - jb * b);

        // Si je suis le propriétaire du bloc (k,jb), je fais la mise à jour fw_row localement.
            if (rank == ownerRow) {
                int localIdx = localIndex[kk * nb + jb];
                if (localIdx != -1) {
                    int* DkJ = &localData[localIdx * blockArea];
                                    // Mise à jour du bloc D(k,J) avec le pivot D(k,k)
                    fw_row(pivotBlock.data(), DkJ, bs, wJ, b);
                                    // Je copie le résultat dans rowBlocks[jb] pour pouvoir le diffuser

                    std::copy(DkJ, DkJ + blockArea, rowBlocks[jb].begin());
                }
            }

              // Ici je lance un broadcast non bloquant du bloc D(k,jb)
        // à partir de son propriétaire ownerRow.
            MPI_Request req;
            MPI_Ibcast(rowBlocks[jb].data(), blockArea, MPI_INT, 
                      ownerRow, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }

    // J'attends que toutes les diffusions de la ligne soient finies
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

       // ===== Phase B.2 : COLONNE de blocs (ib, k) =====
    // Même idée mais pour la colonne : pour chaque bloc (ib,k),
    // le propriétaire fait fw_col puis on Ibcast.
        requests.clear();
        
        for (int ib = 0; ib < nb; ++ib) {
            if (ib == kk) continue; // on saute encore le pivot

            int ownerCol = ownerOf(ib, kk, Pr, Pc);
            int hI = std::min(b, n - ib * b);
        // Si je possède le bloc (ib,k), je fais sa mise à jour fw_col.

            if (rank == ownerCol) {
                int localIdx = localIndex[ib * nb + kk];
                if (localIdx != -1) {
                    int* Dik = &localData[localIdx * blockArea];
                                    // Mise à jour du bloc D(I,k) avec le pivot D(k,k)
                    fw_col(Dik, pivotBlock.data(), hI, bs, b);
                                    // Je copie le résultat dans colBlocks[ib] pour le diffuser

                    std::copy(Dik, Dik + blockArea, colBlocks[ib].begin());
                }
            }
        // Broadcast non bloquant du bloc (ib,k)
            MPI_Request req;
            MPI_Ibcast(colBlocks[ib].data(), blockArea, MPI_INT,
                      ownerCol, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
    // J'attends que toutes les diffusions de la colonne soient finies
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  // ===== Phase C : Blocs internes =====
    // Maintenant que j'ai :
    //  - tous les D(I,k) dans colBlocks[ib]
    //  - tous les D(k,J) dans rowBlocks[jb]
    //
    // Je peux mettre à jour tous les blocs (I,J) qui ne sont ni dans la
    // ligne k ni dans la colonne k, en utilisant fw_inner.        for (int idx = 0; idx < numLocal; ++idx) {
          for (int idx = 0; idx < numLocal; ++idx) { 
    BlockInfo& info = localBlocks[idx];
            int ib = info.bi;
            int jb = info.bj;
 // Si le bloc est dans la ligne ou colonne du pivot,
        // on l'a déjà traité avant.
            if (ib == kk || jb == kk) continue;

            int hI = std::min(b, n - ib * b);
            int wJ = std::min(b, n - jb * b);

            int* Dij = &localData[idx * blockArea];  // bloc (I,J) que je mets à jour
            const int* Dik = colBlocks[ib].data();   // bloc (I,k)
            const int* DkJ = rowBlocks[jb].data();  // bloc (k,J)

                    // Mise à jour complète du bloc interne avec les chemins passant par k
            fw_inner(Dik, DkJ, Dij, hI, wJ, bs, b);
        }
    }

    // ===== Rassemblement =====
    // A la fin, je dois reconstruire la matrice complète n x n sur le rang 0.
// Les autres processus gardent juste leurs blocs locaux dans localData.
    int* D_final = nullptr;
    if (rank == 0) {
        D_final = new int[n * n];
         // Je remplis tout avec INF par défaut, et après je vais écrire
    // les vraies valeurs bloc par bloc.
        for (int i = 0; i < n * n; ++i) D_final[i] = INF;
    }
// Ici je vais passer sur TOUS les rangs r, un par un.
// Pour chaque rang, je récupère la liste des blocs qui lui appartiennent
// (en recalculant computeLocalBlocks avec "r" comme propriétaire).
    for (int r = 0; r < size; ++r) {
        vector<BlockInfo> blocks_r = computeLocalBlocks(n, b, Pr, Pc, r);

        for (size_t idx_r = 0; idx_r < blocks_r.size(); ++idx_r) {
            const BlockInfo& info = blocks_r[idx_r];
            vector<int> buf(blockArea);
            
                    // Si je suis le rang r, je vais copier mon bloc local dans buf.
            if (rank == r) {
                int localIdx = localIndex[info.bi * nb + info.bj];
                if (localIdx != -1) {
                    int* src = &localData[localIdx * blockArea];
                    std::copy(src, src + blockArea, buf.begin());
                }
            }

        // Si r == 0, pas besoin d'envoyer/recevoir, on reste local.
        // Si r != 0, alors le rang r envoie son bloc à 0,
        // et le rang 0 le reçoit dans buf.
            if (r == 0) {
                if (rank == 0) {}
                                // rien à faire de spécial, buf est déjà rempli

            } else {
                if (rank == r) {
                    MPI_Send(buf.data(), blockArea, MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
                if (rank == 0) {
                    MPI_Recv(buf.data(), blockArea, MPI_INT, r, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
// Une fois que buf contient les données du bloc (que ce soit localement
        // ou reçu d'un autre rang), le rang 0 recopie ce bloc dans la vraie
        // matrice D_final, aux bonnes positions globales.
            if (rank == 0) {
                int i0 = info.offset_i;
                int j0 = info.offset_j;
                for (int ii = 0; ii < b; ++ii) {
                    int gi = i0 + ii;   // indice global en ligne
                    if (gi >= n) break;  // on coupe si ça dépasse la vraie matrice
                    for (int jj = 0; jj < b; ++jj) {
                        int gj = j0 + jj; // indice global en colonne
                        if (gj >= n) break; // pareil pour les colonnes
                     // On écrit la valeur du bloc dans la matrice finale

                        D_final[gi * n + gj] = buf[ii * b + jj];
                    }
                }
            }
        }
    }
// Au final, seul le rang 0 renvoie la matrice complète.
// Les autres processus renvoient nullptr, ils n'en ont pas besoin.
    return (rank == 0) ? D_final : nullptr;
}
