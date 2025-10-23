#include "metodoSimplex.h"
#include <stdlib.h>
#include <stdio.h>
#define inf 1.0/0.0     //Infinito, usado para simplificar algumas comparacoes
#define epsilon 1e-12   //Um numero pequeno que e considerado como zero para tornar comparacoes mais seguras

double resolvedorDeProblemaDePL(double** A, double* b, double* c, unsigned long long int m, unsigned long long int n, double* xSaida){
    unsigned long long int* B = (unsigned long long int*) malloc(m*sizeof(unsigned long long int)); //Guarda as colunas da matriz basica
    unsigned long long int* N = (unsigned long long int*) malloc(n*sizeof(unsigned long long int)); //Guarda as colunas da matriz nao basica
    double* lambda = (double*) malloc(m*sizeof(double));    //Guarda o vetor multiplicador Simplex
    double* xB = (double*) malloc(m*sizeof(double));        //Guarda os valores das variaveis na solucao atual
    double* cReduzido = (double*) malloc(n*sizeof(double)); //Guarda os custos reduzidos das colunas nao basicas da solucao atual
    unsigned int exit = 0;          //Flag que controla os casos de saida
    unsigned long long int entra;   //Guarda o indice no vetor que guarda os indices das colunas nao basicas que vai entrar
    unsigned long long int sai;     //Mesma coisa que o anterior mas para o que vai sair da matriz basica
    double min;     //Variavel auxiliar
    unsigned long long int i;       //variavel auxiliar para iteracao
    unsigned long long int j;       //variavel auxiliar para iteracao
    unsigned long long int k;       //variavel auxiliar para iteracao
    unsigned long long int maxCol;  //variavel auxiliar para reoslver sistemas lineares
    double maxVal;  //variavel auxiliar para resolver sistemas lineares
    double val;     //variavel auxiliar par aresolver sistemas lineares

    double** MatrizTemp = (double**) malloc(m*sizeof(double*)); //Matriz que guarda todos os valores de matrizes que tem que ser resolvidas por eliminacao gausiana
    for(i = 0;i < m; ++i){
        MatrizTemp[i] = (double*) malloc(m*sizeof(double));
    }

    double* vetorTemp = (double*) malloc(m*sizeof(double));     //Vetor que guarda todos os valores de vetores uqe tem que serem encontrados por eliminacao gausiana

    for(i = 0; i < n; ++i){ //Inicializa a matriz nao basica com todas as colunas da matriz do problema original
        N[i] = i;
    }

    for(i = 0; i < m; ++i){ //Inicializa a matriz basica com as variaveis ficticias (valores >= n na primeira faze serao sempre as variaveis ficticias)
        B[i] = i + n;
    }

    while(exit == 0){       //Primeira fase
        for(i = 0; i < m; ++i){     //Inicializa a matriz com os valores da matriz basica B
            for(j = 0; j < m; ++j){
                if(B[j] >= n){
                    MatrizTemp[i][j] = (double)(i == (B[j]-n)); //Feito para nao ter que guardar as colunas das variaveis ficticias
                }

                else{
                    MatrizTemp[i][j] = A[i][B[j]];
                }
            }
        }

        for(i = 0; i < m; ++i){     //Inicializa o vetor com os valores de b
            vetorTemp[i] = b[i];
        }

        for(i = 0; i < m; ++i){     //Resolve o sistema por eliminacao gaussiana com pivoteamento para encontrar a solucao basica
            maxCol = i;
            maxVal = (MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i];
            for(j = i + 1; j < m; ++j){
                val = (MatrizTemp[j][i] >= 0) ? MatrizTemp[j][i] : -MatrizTemp[j][i];
                if(val > maxVal){
                    maxVal = val;
                    maxCol = j;
                }
            }

            if(maxCol != i){
                for(j = 0; j < m; ++j){
                    val = MatrizTemp[i][j];
                    MatrizTemp[i][j] = MatrizTemp[maxCol][j];
                    MatrizTemp[maxCol][j] = val;
                }
                val = vetorTemp[i];
                vetorTemp[i] = vetorTemp[maxCol];
                vetorTemp[maxCol] = val;
            }

            if(((MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i]) < epsilon){
                continue;
            }

            val = MatrizTemp[i][i];
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] /= val;
            }

            vetorTemp[i] /= val;
            for(j = 0; j < m; ++j){
                if(j == i){
                    continue;
                }

                val = MatrizTemp[j][i];
                for(k = 0; k < m; ++k){
                    MatrizTemp[j][k] -= val * MatrizTemp[i][k];
                }

                vetorTemp[j] -= val * vetorTemp[i];
            }
        }

        for(i = 0; i < m; ++i){     //Copia do vetor temporario para o vetor que guarda a solucao basica
            xB[i] = vetorTemp[i];
        }

        for(i = 0; i < m; ++i){     //Inicializa a matriz com B transposto
            for(j = 0; j < m; ++j){
                if(B[j] >= n){
                    MatrizTemp[j][i] = (double)(i == (B[j]-n)); //Feito para nao ter que guardar as colunas das variaveis ficticias
                }

                else{
                    MatrizTemp[j][i] = A[i][B[j]];
                }
            }
        }

        for(i = 0; i < m; ++i){     //Inicializa o vetor com os custos basicos (para a primeira faze, sao todos 0 exceto pelos ficticios)
            if(B[i] >= n){
                vetorTemp[i] = 1;
            }

            else{
                vetorTemp[i] = 0;
            }
        }

        for(i = 0; i < m; ++i){     //Resolve o sistema por eliminacao gaussiana com pivoteamento para encontrar o vetor multiplicador simplex
            maxCol = i;
            maxVal = (MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i];
            for(j = i + 1; j < m; ++j){
                val = (MatrizTemp[j][i] >= 0) ? MatrizTemp[j][i] : -MatrizTemp[j][i];
                if(val > maxVal){
                    maxVal = val;
                    maxCol = j;
                }
            }
            if(maxCol != i){
                for(j = 0; j < m; ++j){
                    val = MatrizTemp[i][j];
                    MatrizTemp[i][j] = MatrizTemp[maxCol][j];
                    MatrizTemp[maxCol][j] = val;
                }

                val = vetorTemp[i];
                vetorTemp[i] = vetorTemp[maxCol];
                vetorTemp[maxCol] = val;
            }

            if(((MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i]) < epsilon){
                continue;
            }

            val = MatrizTemp[i][i];
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] /= val;
            }

            vetorTemp[i] /= val;
            for(j = 0; j < m; ++j){
                if(j == i){
                    continue;
                }

                val = MatrizTemp[j][i];
                for(k = 0; k < m; ++k){
                    MatrizTemp[j][k] -= val * MatrizTemp[i][k];
                }

                vetorTemp[j] -= val * vetorTemp[i];
            }
        }

        for(i = 0; i < m; ++i){     //Copia do vetor temporario para o vetor multiplicador simplex
            lambda[i] = vetorTemp[i];
        }

        min = inf;                  //Nada pode ser maior que infinito :P (feito por preguica para nao ter que usar o primeiro valor)

        for(i = 0; i < n; ++i){     //Calcula os custos reduzidos das colunas nao basicas (e encontra a coluna que vai entrar na base)
            if(N[i] < n){
                val = 0;
                for(j = 0; j < m; ++j){
                    if(A[j][N[i]] != 0 && lambda[j] != 0){
                        val += lambda[j] * A[j][N[i]];
                    }
                }
                cReduzido[i] = -val;
            }

            else{
                cReduzido[i] = 1 - lambda[N[i] - n];
            }

            if(min >= cReduzido[i]){
                entra = i;
                min = cReduzido[i];
            }
        }

        if(min >= 0){               //Se o menor for positivo, todos sao positivos e estamos no otimo
            exit = 1;
            break;
        }

        for(i = 0; i < m; ++i){     //Inicializa a matriz com B
            for(j = 0; j < m; ++j){
                if(B[j] >= n){
                    MatrizTemp[i][j] = (double)(i == (B[j]-n)); //Feito para nao ter que guardar as colunas das variaveis ficticias
                }

                else{
                    MatrizTemp[i][j] = A[i][B[j]];
                }
            }
        }

        if(N[entra] >= n){          //Inicializa com a coluna nao basica que vai entrar (caso de ser variavel ficticia)
            for(i = 0; i < m; ++i){
                vetorTemp[i] = (double)(i == (N[entra]-n)); //Feito para nao ter que guardar as colunas das variaveis ficticias
            }
        }

        else{                       //Inicializa com a coluna nao basica que vai entrar (caso de ser variavel nao ficticia)
            for(i = 0; i < m; ++i){
                vetorTemp[i] = A[i][N[entra]];
            }
        }

        for(i = 0; i < m; ++i){     //Resolve o sistema por eliminacao gaussiana com pivoteamento para achar os y
            maxCol = i;
            maxVal = (MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i];
            for(j = i + 1; j < m; ++j){
                val = (MatrizTemp[j][i] >= 0) ? MatrizTemp[j][i] : -MatrizTemp[j][i];
                if(val > maxVal){
                    maxVal = val;
                    maxCol = j;
                }
            }
            if(maxCol != i){
                for(j = 0; j < m; ++j){
                    val = MatrizTemp[i][j];
                    MatrizTemp[i][j] = MatrizTemp[maxCol][j];
                    MatrizTemp[maxCol][j] = val;
                }

                val = vetorTemp[i];
                vetorTemp[i] = vetorTemp[maxCol];
                vetorTemp[maxCol] = val;
            }

            if(((MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i]) < epsilon){
                continue;
            }

            val = MatrizTemp[i][i];
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] /= val;
            }

            vetorTemp[i] /= val;
            for(j = 0; j < m; ++j){
                if(j == i){
                    continue;
                }

                val = MatrizTemp[j][i];
                for(k = 0; k < m; ++k){
                    MatrizTemp[j][k] -= val * MatrizTemp[i][k];
                }

                vetorTemp[j] -= val * vetorTemp[i];
            }
        }

        min = inf;                  //De novo, nada pode ser maior que infinito
        exit = 2;
        for(i = 0; i < m; ++i){     //Calcula o tamanho de passo e confere se o problema e ilimitado (descobre a coluna que sai da base)
            if(vetorTemp[i] > 0){
                exit = 0;
                val = xB[i]/vetorTemp[i];
                if(min > val){
                    sai = i;
                    min = val;
                }
            }
        }

        if(exit == 2){              //Se isso for o caso, quer dizer que todos sao negativos e o problema e ilimitado
            break;
        }

                                    //A coluna que sai, sai e a que entra, entra
        i = B[sai];
        B[sai] = N[entra];
        N[entra] = i;
    }

    for(i = 0; i < m; ++i){ //Remove as variaveis artificiais de B (se achar alguma na base igual a 0 muda para alguma das não basicas para dar certo, se não for 0 é porque é infactivel)
        if(B[i] >= n){      //Verifica se B[i] é variável artificial
            if(((xB[i] < 0) ? -xB[i] : xB[i]) > epsilon){ //Alguma xB está na base e é não zero no ótimo (problema impossivel)
                exit = 2;
                break;
            }

            for(j = 0; j < n; ++j){ // Substitui por qualquer variável não básica com valor zero
                if(N[j] < n){
                    maxCol = B[i];
                    B[i] = N[j];
                    N[j] = maxCol;
                    break;
                }
            }
        }
    }

    if(exit == 2){          //Problema Impossível
        for(i = 0; i < n; ++i){
            xSaida[i] = 0.0/0.0;
        }
        min = 0.0/0.0;
        exit = 3;
    }

    else{                   //Problema nao impssivel, e a matriz basica e uma solucao basica factivel, pode comecar a segunda faze
        exit = 0;
    }

    if(exit == 0){          //Remove as variaveis artificiais da matriz não basica (teoricamente a unica que deve ter artificiais depois do passo de tirar elas da basica... eu espero :p O que isso já deu de problema é brincadeira) e printa a base inicial factivel
        j = 0;
        for(i = 0; i < n; ++i){
            if(N[i] < n){
                N[j] = N[i];
                j++;
            }
        }
                            //Se for realmente querer usar esse programa para alguma aplicação mais de verdade, remover essa parte em particular
        printf("\nPartição basica inicial factível:\nB = [");
        for(i = 0; i < m; ++i){
            printf("%llu",B[i] + 1);
            if(i < m-1){
                printf(", ");
            }
        }
        printf("]");
    }

    while(exit == 0){       //Segunda Fase. B agora guarda uma particao basica factível
        for(i = 0; i < m; ++i){     //Inicializa a matriz com os valores da matriz basica B
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] = A[i][B[j]];
            }
        }

        for(i = 0; i < m; ++i){     //inicializa o vetor com os valores de b
            vetorTemp[i] = b[i];
        }

        for(i = 0; i < m; ++i){     //Resolve o sistema por eliminacao gaussiana com pivoteamento para encontrar a solucao basica
            maxCol = i;
            maxVal = (MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i];
            for(j = i + 1; j < m; ++j){
                val = (MatrizTemp[j][i] >= 0) ? MatrizTemp[j][i] : -MatrizTemp[j][i];
                if(val > maxVal){
                    maxVal = val;
                    maxCol = j;
                }
            }

            if(maxCol != i){
                for(j = 0; j < m; ++j){
                    val = MatrizTemp[i][j];
                    MatrizTemp[i][j] = MatrizTemp[maxCol][j];
                    MatrizTemp[maxCol][j] = val;
                }

                val = vetorTemp[i];
                vetorTemp[i] = vetorTemp[maxCol];
                vetorTemp[maxCol] = val;
            }

            if(((MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i]) < epsilon){
                continue;
            }

            val = MatrizTemp[i][i];
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] /= val;
            }

            vetorTemp[i] /= val;
            for(j = 0; j < m; ++j){
                if(j == i){
                    continue;
                }

                val = MatrizTemp[j][i];
                for(k = 0; k < m; ++k){
                    MatrizTemp[j][k] -= val * MatrizTemp[i][k];
                }

                vetorTemp[j] -= val * vetorTemp[i];
            }
        }

        for(i = 0; i < m; ++i){     //Copia do vetor temporario para o vetor que guarda a solucao basica
            xB[i] = vetorTemp[i];
        }

        for(i = 0; i < m; ++i){     //Inicializa a matriz com B transposto
            for(j = 0; j < m; ++j){
                MatrizTemp[j][i] = A[i][B[j]];
            }
        }

        for(i = 0; i < m; ++i){     //Inicializa o vetor com os custos basicos
            vetorTemp[i] = c[B[i]];
        }

        for(i = 0; i < m; ++i){     //Resolve o sistema por eliminacao gaussiana com pivoteamento para encontrar o vetor multiplicador simplex
            maxCol = i;
            maxVal = (MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i];
            for(j = i + 1; j < m; ++j){
                val = (MatrizTemp[j][i] >= 0) ? MatrizTemp[j][i] : -MatrizTemp[j][i];
                if(val > maxVal){
                    maxVal = val;
                    maxCol = j;
                }
            }

            if(maxCol != i){
                for(j = 0; j < m; ++j){
                    val = MatrizTemp[i][j];
                    MatrizTemp[i][j] = MatrizTemp[maxCol][j];
                    MatrizTemp[maxCol][j] = val;
                }

                val = vetorTemp[i];
                vetorTemp[i] = vetorTemp[maxCol];
                vetorTemp[maxCol] = val;
            }

            if(((MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i]) < epsilon){
                continue;
            }

            val = MatrizTemp[i][i];
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] /= val;
            }

            vetorTemp[i] /= val;
            for(j = 0; j < m; ++j){
                if(j == i){
                    continue;
                }

                val = MatrizTemp[j][i];
                for(k = 0; k < m; ++k){
                    MatrizTemp[j][k] -= val * MatrizTemp[i][k];
                }

                vetorTemp[j] -= val * vetorTemp[i];
            }
        }

        for(i = 0; i < m; ++i){     //Copia do vetor temporario para o vetor multiplicador simplex
            lambda[i] = vetorTemp[i];
        }

        min = inf;                  //Inicializa com infinito por preguica

        for(i = 0; i < n-m; ++i){   //Calcula os custos reduzidos das colunas nao basicas (e encontra a coluna que vai entrar na base)
            val = 0;
            for(j = 0; j < m; ++j){
                if(A[j][N[i]] != 0 && lambda[j] != 0){
                    val += lambda[j] * A[j][N[i]];
                }
            }

            cReduzido[i] = c[N[i]] - val;
            if(min > cReduzido[i]){
                entra = i;
                min = cReduzido[i];
            }
        }

        if(min >= 0){               //Se o menor for positivo, todos sao positivos e estamos no otimo
            exit = 1;
            break;
        }

        for(i = 0; i < m; ++i){     //Inicializa a matriz com B
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] = A[i][B[j]];
            }
        }

        for(i = 0; i < m; ++i){     //Inicializa com a coluna nao basica que vai entrar
            vetorTemp[i] = A[i][N[entra]];
        }

        for(i = 0; i < m; ++i){     //Resolve o sistema por eliminacao gaussiana com pivoteamento para achar os y
            maxCol = i;
            maxVal = (MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i];
            for(j = i + 1; j < m; ++j){
                val = (MatrizTemp[j][i] >= 0) ? MatrizTemp[j][i] : -MatrizTemp[j][i];
                if(val > maxVal){
                    maxVal = val;
                    maxCol = j;
                }
            }

            if(maxCol != i){
                for(j = 0; j < m; ++j){
                    val = MatrizTemp[i][j];
                    MatrizTemp[i][j] = MatrizTemp[maxCol][j];
                    MatrizTemp[maxCol][j] = val;
                }

                val = vetorTemp[i];
                vetorTemp[i] = vetorTemp[maxCol];
                vetorTemp[maxCol] = val;
            }

            if(((MatrizTemp[i][i] >= 0) ? MatrizTemp[i][i] : -MatrizTemp[i][i]) < epsilon){
                continue;
            }

            val = MatrizTemp[i][i];
            for(j = 0; j < m; ++j){
                MatrizTemp[i][j] /= val;
            }

            vetorTemp[i] /= val;
            for(j = 0; j < m; ++j){
                if(j == i){
                    continue;
                }

                val = MatrizTemp[j][i];
                for(k = 0; k < m; ++k){
                    MatrizTemp[j][k] -= val * MatrizTemp[i][k];
                }

                vetorTemp[j] -= val * vetorTemp[i];
            }
        }

        min = inf;                  //Preguica
        exit = 2;
        for(i = 0; i < m; ++i){     //Calcula o tamanho de passo e confere se o problema e ilimitado (descobre a coluna que sai da base)
            if(vetorTemp[i] > 0){
                exit = 0;
                val = xB[i]/vetorTemp[i];
                if(min > val){
                    sai = i;
                    min = val;
                }
            }
        }

        if(exit == 2){              //Se isso for o caso, quer dizer que todos sao negativos e o problema e ilimitado
            break;
        }

        printf("\n\nPartição basica atual:\n");
        for(i = 0; i < m; ++i){
            printf("x%llu = %lf ",B[i] + 1,xB[i]);
        }
        
                                    //A coluna que sai, sai e a que entra, entra
        i = B[sai];
        B[sai] = N[entra];
        N[entra] = i;
    }

    if(exit == 1){          //Se exit for 1, indica que o problema tem solucao otima
        min = 0;
        for(i = 0; i < m; ++i){     //Calcula o valor otimo
            min = min + xB[i]*c[B[i]];
        }

        for(i = 0; i < n ; ++i){    //Escreve no vetor de saida
            xSaida[i] = 0;
        }
        for(j = 0; j < m; ++j){
            if(B[j] >= n){
                continue;
            }
            xSaida[B[j]] = xB[j];
        }
    }

    else if(exit == 2){    //Se exit for 2, indica que o problema e ilimitado
        min = -inf;
        for(i = 0; i < n; ++i){
            xSaida[i] = 0.0/0.0;
        }
    }

    free(B);
    free(N);
    free(lambda);
    free(xB);
    free(cReduzido);
    for(i = 0; i < m; ++i){
        free(MatrizTemp[i]);
    }

    free(MatrizTemp);
    free(vetorTemp);

    return min;
}

#undef inf
#undef epsilon

#include <stdio.h>
#include <stdlib.h>

void promptResolvedorDePL(){
    unsigned long long int n0; 
    unsigned long long int m0;
    unsigned long long int i;
    unsigned long long int j;

    printf("Número de variáveis: ");
    if(scanf("%llu", &n0) != 1 || n0 == 0){
        printf("Entrada inválida.\n");
        return;
    }

    printf("Número de restrições: ");
    if(scanf("%llu", &m0) != 1 || m0 == 0){
        printf("Entrada inválida.\n");
        return;
    }

    double **A = (double**) malloc(m0 * sizeof(double*));
    for(i = 0; i < m0; i++){
        A[i] = (double*) malloc(n0 * sizeof(double));
        for(j = 0; j < n0; j++){
             A[i][j] = 0.0;
        }
    }

    double* b = (double*) malloc(m0 * sizeof(double));
    double* c = (double*) malloc(n0 * sizeof(double));
    char** rels = (char**) malloc(m0 * sizeof(char*));
    for(i = 0; i < m0; i++){ 
        rels[i] = (char*)malloc(3 * sizeof(char));
    }

    for(i = 0; i < m0; i++){
        b[i] = 0.0;
    }

    for(j = 0; j < n0; j++){
        c[j] = 0.0;
    }

    char tipo_objetivo[8];
    printf("\nFunção objetivo (digite 'min' ou 'max'): ");
    scanf("%7s", tipo_objetivo);
    int maximizar = 0;
    if(tipo_objetivo[0] == 'M' || tipo_objetivo[0] == 'm'){
        if(tipo_objetivo[1] == 'A' || tipo_objetivo[1] == 'a'){
            maximizar = 1;
        }
    }

    printf("\nCoeficientes da função objetivo:\n");
    for(j = 0; j < n0; j++){
        printf("c%llu = ", j + 1);
        if(scanf("%lf", &c[j]) != 1){
            c[j] = 0.0;
        } 
    }

    if(maximizar){
        for(j = 0; j < n0; j++){
            c[j] = -c[j];
        }
    }

    printf("\nRestrições:\n");
    for(i = 0; i < m0; i++){
        for(j = 0; j < n0; j++){
            printf("A(%llu,%llu) = ", i + 1, j + 1);
            if(scanf("%lf", &A[i][j]) != 1){ 
                A[i][j] = 0.0;
            }
        }

        printf("Relação (<=, >=, =): ");
        scanf("%2s", rels[i]);
        printf("b%llu = ", i + 1);
        if(scanf("%lf", &b[i]) != 1){
            b[i] = 0.0;
        }
    }

    unsigned long long folgas = 0; //Muda para a forma padrão
    for(i = 0; i < m0; i++){
        if(b[i] < 0.0){ //Gabrante b >= 0
            for(j = 0; j < n0; j++){ 
                A[i][j] = -A[i][j];
            }

            b[i] = -b[i];
            if(rels[i][0] == '<' && rels[i][1] == '='){
                rels[i][0] = '>'; 
                rels[i][1] = '=';
            }

            else if(rels[i][0] == '>' && rels[i][1] == '='){
                rels[i][0] = '<'; 
                rels[i][1] = '='; 
            }
        }

        if((rels[i][0] == '<' && rels[i][1] == '=') || (rels[i][0] == '>' && rels[i][1] == '=')){   
            folgas++;
        }
    }

    unsigned long long n_total = n0 + folgas;
    double** A_std = (double**) malloc(m0 * sizeof(double*)); 
    for(i = 0; i < m0; i++){
        A_std[i] = (double*) malloc(n_total * sizeof(double));
        for(j = 0; j < n_total; j++){
            A_std[i][j] = 0.0;
        }
    }
 
    unsigned long long col = n0; //Bota o A na forma padrão e guarda no A_std
    for(i = 0; i < m0; i++){
        for(j = 0; j < n0; j++){
            A_std[i][j] = A[i][j];
        }

        if(rels[i][0] == '<' && rels[i][1] == '='){
            A_std[i][col++] = 1.0;
        } 
        
        else if(rels[i][0] == '>' && rels[i][1] == '='){
            A_std[i][col++] = -1.0;
        }
    }

    double* c_ext = (double*) malloc(n_total * sizeof(double));
    for(j = 0; j < n_total; j++){
        if(j < n0){ 
            c_ext[j] = c[j];
        }

        else{ 
            c_ext[j] = 0.0;
        }
    }

    printf("\nProblema em forma padrão:\n");
    printf("Min z = ");
    for(j = 0; j < n0; j++){
        if(j > 0){
            if(c[j] >= 0){ 
                printf(" + ");
            }
            
            else{
                printf(" - ");
            }
        }

        printf("%.6g x%llu", c[j] >= 0 ? c[j] : -c[j], j + 1);
    }

    printf("\n\nSujeito a:\n");
    int primeiro;
    for(i = 0; i < m0; i++){
        primeiro = 1;
        for(j = 0; j < n_total; j++){
            if(A_std[i][j] != 0.0){
                if(!primeiro){
                    if(A_std[i][j] >= 0){
                        printf(" + ");
                    }

                    else{
                        printf(" - ");
                    }
                } 

                else if(A_std[i][j] < 0){
                    printf("-");
                }

                printf("%.6g x%llu", A_std[i][j] >= 0 ? A_std[i][j] : -A_std[i][j], j + 1);
                primeiro = 0;
            }
        }

        printf(" = %.6g\n", b[i]);
    }

    printf("\nTodas as variáveis xj >= 0.\n");
    double *xsaida = (double *)malloc(n_total * sizeof(double)); //inicializa o vetor de respostas
    for(j = 0; j < n_total; j++){
        xsaida[j] = 0.0;
    }

    double valor_obj = resolvedorDeProblemaDePL(A_std, b, c_ext, m0, n_total, xsaida); //Finalmente resolver o problema de PL :D
    if(!(valor_obj == valor_obj)){
        printf("\nProblema inviável ou sem solução.\n");
    }
    
    else if(valor_obj > 1e200 || valor_obj < -1e200){
        printf("\nProblema ilimitado.\n");
    } 
    
    else{
        printf("\n\nValor ótimo (min) = %.12g\n", valor_obj);
        printf("Solução:\n");
        for(j = 0; j < n0; j++){
            printf("x%llu = %.12g\n", j + 1, xsaida[j]);
        }
    }

    free(xsaida);
    free(c_ext);
    for(i = 0; i < m0; i++){
        free(A_std[i]);
    }

    free(A_std);
    for(i = 0; i < m0; i++){
        free(A[i]);
    }

    free(A);
    for(i = 0; i < m0; i++){
        free(rels[i]);
    }
    free(rels);
    free(b);
    free(c);
}
