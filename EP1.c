#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define eps 10e-5

//////////////////////////////////////////////////////////////////////
///////														                             ///////
/////// Rodrigo Fill Rangel - NUSP : 4359840     		       	   ///////
///////														                             ///////
/////// Roberto Costa Ceccato - NUSP : 4103944			           ///////
///////												                               	 ///////
/////// Turma 03 - Linguagem: C                                ///////
//////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------FUNÇÕES DE MATRIZES-----------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

// Função responsável por alocar dinâmicamente as matrizes usadas pelo programa, definimos estas como um ponteiro //
// duplo do tipo double, para os cálculos futuros. //

double** CriarMatriz(int n, int m) {
    double** Matriz = malloc (n * sizeof (double*));

    for (int i = 0; i < n; i++ ){
        Matriz[i] = (double*)malloc(m * sizeof(double));
        for (int j = 0; j < m; j++)
          Matriz[i][j] = 0;
    }

    
 return Matriz;
}

// DeletarMatriz é a função que usa o free da bíblioteca stdlib, de forma a, como o nome diz, deletar as matrizes  //
// criadas que não serão mais usadas. Esta função é uma importante forma de gerenciar o uso de memória RAM do compu-//
// tador, que é bastante importante neste programa. //

void DeletarMatriz (double** M, double n) {
    int i;

    for(i = 0; i < n; i++)
      free(M[i]);
    free(M);
}

// Esta função é responsável por imprimir na tela os valores armazenados no ponteiro duplo passado como parâmetro, //
// usamos três digitos depois da vírgula, assumindo essa como a significância dos problemas, além de garantir que as //
// matrizes não fiquem desnecessáriamente grandes em tela. //

void printMatriz (double** Z, int n, int m) {

    for (int e = 0; e < 6*m + 3; e++)
      printf("-");
    printf ("\n");

    for ( int i = 0; i < n; i++ ) {
      for ( int j = 0; j < m; j++ ) {
        if(j == 0 && m == 1)
          printf ("| %.3f |", Z[i][j]);
        else if (j == 0)
          printf ("| %.3f", Z[i][j]);
        else if (j == m - 1)
          printf (" %.3f |", Z[i][j]);
        else
          printf (" %.3f", Z[i][j]);
      }
      printf ("\n");
    }

    for (int e = 0; e < 6*m + 3; e++)
        printf("-");
    printf ("\n");
}

// Função responsável por transpor as linhas e colunas da matriz W passada como parâmetro, esta função cria uma matriz //
// T, que será retornada como a transposta de W. //

void TransporMatriz (double** de, double** para, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            para[j][i] = de[i][j];
        }
    }
}

// A função abaixo copia a matriz A passada como parâmetro na matriz B, necessário nas tarefas subsequentes. //

void Copiar_Matriz(double** B, double** A, int n,int m){
    for (int i = 0; i < n; i++){
      for (int j = 0; j < m; j++)
        B[i][j] = A[i][j];
    }
}

// Esta função é responsável por inicializar as matrizes criadas, de forma a evitar problemas com lixo armazenado //
// nas posições dos ponteiros alocados, ela atribui zero a todas as posições da matriz. //

void Inicializar_Zeros(double** W, int n, int m){
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++) {
      W[i][j] = 0;
    }
  }
}

// Abaixo está a função que fará a multiplicação de duas matrizes, passadas como parâmetro, ela retorna a matriz A, //
// produto das matrizes W e H. //

void Multiplicar_Matriz (double** W, double** H, double ** res, int n, int m, int p) {
  int i, j, k;
  double soma;

  for (i = 0; i < n; i++ ){
    for (j = 0; j < m; j++ ){
      soma = 0;
      for (k = 0; k < p; k++ ){
      soma += W[i][k]*H[k][j];
      }
      res[i][j] = soma;
    }
  }
}

//----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------PRIMEIRA TAREFA--------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

// A primeira função da primeira tarefa fará o calculo do seno e cosseno da rotação de givens. Para tal serão passa-//
// dos dois parâmetros, Wik e Wjk, que serão os valores armazenados nas respectivas posições das matrizes às quais  //
// será aplicado a rotação de givens. Os outros dois parâmetros, c e s, são ponteiros usados como forma de retornar //
// ambos os valores de seno e cosseno. //
void Calc_Sen_Cos_QR(double Wik, double Wjk, double *c, double *s) {
    double T;                        // Variável auxiliar para cálculo do seno e do cosseno. //

    if ( fabs(Wik) > fabs(Wjk) ) {   // Condição para cálculo de seno e cosseno da forma descrita abaixo. //
        T = - Wjk/Wik;
        *c = 1/sqrt(1 + T*T);
        *s = *c * T;
    }
    else {                           // Outra forma para cálculo de seno e cosseno. //
        T = - Wik/Wjk;
        *s = 1/sqrt(1 + T*T);
        *c = *s * T;
    }
}

// Função que calcula a rotação de givens à uma linha da matriz passada como parâmetro. Como esta é uma das principais//
// funções do exercício programa, pois durante grande parte do tempo o programa estará executando esta função, //
// a implementação é identica à especificada no enúnciado //
void Rotacao_Givens ( double** W, int n, int m, int i, int j, double c, double s) {
  double temp;

      for ( int k = 0; k < m; k++ ) {
      temp = c * W[i][k] - s * W[j][k];
      W[j][k] = s * W[i][k] + c * W[j][k];
      W[i][k] = temp;
     }
}

// A terceira função da primeira tarefa também é implementada de forma identica à presente no enunciado do EP, esta //
// função deve calcular a matriz X referente à resposta da equação W*H = b, por meio de várias rotações de givens,  //
// ao longo de toda a matriz W e do vetor b, tornando W uma matriz diagonal superior. //

void Decomposicao_QR(double** W, int n, int m, double** X, double** b) {

    double somatoria, c, s;
    int i, j, k;
    for ( k = 0; k < m ; k++ ) {
        for ( j = n-1; j > k ; j-- ) {
            i = j - 1;                                      // A rotação de givens depende da linha anterior //
            if ( fabs(W[j][k]) > eps ) {
                Calc_Sen_Cos_QR(W[i][k], W[j][k], &c, &s);  // Calculo de sen e cos para W[i][k] e W[j][k] //
                Rotacao_Givens(W, n, m, i, j, c, s);        // Aplica Q(i,j,θ) à matriz W (com c e s definidos acima) // 
                Rotacao_Givens(b, n, 1, i, j, c, s);        // Aplica Q(i,j,θ) ao vetor b (com c e s definidos acima) // 
            }
        }
    }

    for ( k = m - 1; k >= 0; k-- ) {
      somatoria = 0;                      
      for ( j = k; j < m; j++ )
        somatoria += W[k][j] * X[j][0];                    
      X[k][0] = (b[k][0] - somatoria )/W[k][k];
    }
}

// A última função da primeira tarefa deste EP corresponde à uma implementação análoga à anterior, porém, esta  //
// corresponde, não à um vetor resposta b, mas à uma matriz resposta A, desta forma, faremos o mesmo da função  //
// anterior, porém aplicando a todas as colunas de "b". //

void Decomposicao_QR_Multipla (double** W, int n, int p, int m, double** H, double** A) {
    double sum, c, s;

    for ( int k = 0; k < p; k++ ) {
        for ( int j = n - 1; j > k; j-- ) {
            int i = j - 1;
            if (W[j][k] != 0 ) {
              Calc_Sen_Cos_QR(W[i][k], W[j][k], &c, &s);
              Rotacao_Givens(W, n, p, i, j, c, s);
              Rotacao_Givens(A, n, m, i, j, c, s);
            }
        }
    }
    for ( int k = p - 1; k >= 0; k--  ) {   
    for ( int j = 0; j < m; j++  ) {
      sum = 0;
      for (int i = k ; i < p ; i++){        // Para o cálculo da matriz resposta H teremos mais um for //
        sum += W[k][i] * H[i][j];
      }
      H[k][j] = (A[k][j] - sum)/W[k][k];
    }
  }
}

void Primeira_Tarefa (){
//----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------PRIMEIRA TAREFA a)------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

  double** W = CriarMatriz(64, 64);
  double** a = CriarMatriz(64, 1);
  double** X = CriarMatriz(64, 1);
  Inicializar_Zeros(X, 64, 1);

  // Atribuíndo os valores à matriz W//
  for (int i = 0; i < 64; i ++) {
    for (int j = 0; j < 64; j++ ) {
      if (abs(i-j) == 1)
          W[i][j] = 1;

      if (i == j)
          W[i][j] = 2;

      if ( abs(i - j) > 1)
        W[i][j] = 0;
    }
  }

  // Atribuíndo os valores à matriz b//
  for (int i = 0; i < 64; i++){
    a[i][0] = 1;
  }

  printf("\n");

  // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e < 50; e++){
    if (e == 25)
        printf("PRIMEIRA TAREFA a)");
    else
        printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  printf("\n");
  // fim do cabeçalho //

  Decomposicao_QR(W, 64, 64, X, a);
  printf("Matriz X: \n");
  printMatriz(X, 64, 1);

  printf("\n");

  DeletarMatriz(W, 64);
  DeletarMatriz(a, 64);
  DeletarMatriz(X, 64);

  printf("\n");
  printf("\n");

//----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------PRIMEIRA TAREFA b)------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

  double** M = CriarMatriz(20, 17);
  double** b = CriarMatriz(20, 1);
  double** Y = CriarMatriz(17, 1);
  Inicializar_Zeros(Y, 17, 1);

  // Atribuíndo os valores à matriz M //
  for (int i = 0; i < 20;i ++) {
    for (int j = 0; j < 17; j++ ) {
      if(abs(i-j) <= 4)
        M[i][j] = 1/((double)i+(double)j+1);
      else
        M[i][j] = 0;
    }
  }

  // Atribuíndo os valores à matriz b//
  for (int i = 0; i < 20; i++){
      b[i][0] = 1;
  }

  // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e < 50; e++){
    if (e == 25)
        printf("PRIMEIRA TAREFA b)");
    else
        printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  // Fim do cabeçalho //

  printf("\n");

  Decomposicao_QR(M, 20, 17, Y, b);

  printf("\n");
  printf("Matriz X: \n");
  printMatriz(Y, 17, 1);

  printf("\n");

  DeletarMatriz(M, 20);
  DeletarMatriz(Y, 17);
  DeletarMatriz(b, 20);    

  printf("\n");
  printf("\n");


//----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------PRIMEIRA TAREFA c)------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

  double** V = CriarMatriz(64, 64);       // Matriz fornecida V //
  double** c = CriarMatriz(64, 3);        // Matriz resposta c  //
  double** Z = CriarMatriz(64, 3);        // Matriz incógnita Z //
  Inicializar_Zeros(Z, 64, 3);

  // Atribuíndo os valores à matriz W//
  for (int i = 0; i < 64; i ++) {
    for (int j = 0; j < 64; j++ ) {
      if (abs(i-j) == 1)
        V[i][j] = 1;

      if (i == j)
        V[i][j] = 2;

      if ( abs(i - j) > 1)
        V[i][j] = 0;
    }
  }

  // Atribuíndo os valores à matriz b//
  for (int i = 0; i < 64; i++){
    c[i][0] = 1;
    c[i][1] = i
    ;
    c[i][2] = 2*i - 1;
  }

  printf("\n");

  // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e < 50; e++){
    if (e == 25)
        printf("PRIMEIRA TAREFA c)");
    else
        printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  printf("\n");
  // fim do cabeçalho //

  Decomposicao_QR_Multipla(V, 64, 64, 3, Z, c);

  printf("Matriz H: \n");
  printMatriz(Z, 64, 3);

  printf("\n");

  DeletarMatriz(V, 64);
  DeletarMatriz(c, 64);
  DeletarMatriz(Z, 64);

  printf("\n");
  printf("\n");

//----------------------------------------------------------------------------------------------------------------//
//----------------------------------------------PRIMEIRA TAREFA d)------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

  double** N = CriarMatriz(20, 17); // Matriz fornecida N //
  double** d = CriarMatriz(20, 3);  // Matriz resposta d  //
  double** Q = CriarMatriz(17, 3);  // Matriz incognita Q //
  Inicializar_Zeros(Q, 17, 3);

  for (int i = 0; i < 20;i ++) {
    for (int j = 0; j < 17; j++ ) {
      if(abs(i-j) <= 4){
        N[i][j] = 1/((double)i+(double)j+1);
      }
      else{
        N[i][j] = 0;
      }
    }
  }

  for (int i = 0; i < 20; i++){
    d[i][0] = 1;
    d[i][1] = i + 1;
    d[i][2] = 2*i + 1;
  }

  // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e < 50; e++){
    if (e == 25)
        printf("PRIMEIRA TAREFA d)");
    else
        printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  printf("\n");
  // fim do cabeçalho //

  Decomposicao_QR_Multipla(N, 20, 17, 3, Q, d);

  printf("Matriz H: \n");
  printMatriz(Q, 17, 3);

  printf("\n");

  DeletarMatriz(N, 20);
  DeletarMatriz(d, 20);
  DeletarMatriz(Q, 17);

  printf("\n");
  printf("\n");
}

//----------------------------------------------------------------------------------------------------------------//
//------------------------------------------------SEGUNDA TAREFA--------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

// Para a primeira função da segunda tarefa deste exercicio programa temos a função de normalizar colunas de uma  //
// matriz fornecida como parâmetro, aqui introduzimos uma condição para realizar a normalização, por meio de um   //
// if, que verifica se a soma é nula, caso seja, não realizamos a divisão. Isso porque se a somatoria for nula    //
// teremos a divisão por um zero, resultando em um nan. De resto a implementação é semelhante ao enúnciado. // 
void NormalizarColunas (double** W, int n, int p) {
  double somatoria;

  for ( int j = 0 ; j < p; j++ ) {
  somatoria = 0;
    for ( int i = 0; i < n; i++ ) {
      somatoria += (W[i][j]*W[i][j]);
    }
    if(somatoria > eps){                  // Aqui a condição inserida pela dupla para evitar os nan's //
      for ( int k = 0; k < n; k++ ) {
        W[k][j] /= sqrt(somatoria);
      }
    }  
  }
}

// Para a função abaixo temos a implementação da expressão matemática: max{0, H[i][j]}. Ela é muito relevante para //
// esta parte da programa pois é nesta segunda tarefa que buscamos implementar a fatoração não negativas de matrizes //
// esta função então acaba com os elementos negativos das matrizes. //
void Maximo_0 (double** H, int p, int m) {
  for(int i = 0; i < p; i++){
    for (int j = 0; j < m; j++){
      if (H[i][j] < eps)
        H[i][j] = 0;
    }
  }
}

// A terceira função da segunda tarefa inicializa a matriz Q, passada como parâmetro, de forma aleatória. //
void Randomizador (double ** Q, int n, int m) {
  time_t t;
  srand((unsigned) time(&t));  
  for ( int i = 0; i < n; i++ ) {
      for ( int j = 0; j < m; j++ ) {
          Q[i][j] = (double)rand();
      }
  }
}

// Esta função é responsável por fazer a fatoração não negativa de uma matriz, conforme especificado pelo item 3 do  //
// enunciado do EP. Ela recebe como parâmetro uma matriz A, aquela que dará origem às outras duas matrizes fatoradas,// 
// W e H, recebe também os valores n, m e p, sendo que o primeiro corresponde ao número de linhas das matrizes A e W,//
// o segundo corresponde ao número de colunas da matriz W, mas também corresponde ao número de linhas da matriz H.   //
// Por fim, o último corresponde ao número de colunas de W e de H. A seguir comentaremos melhor cada passagem.//
double** Fatoracao_nao_Negativa(double** A, int n, int p, int m){
    double** W = CriarMatriz(n, p);        // Cria a matriz W, que será retornada ao final//
    Randomizador(W, n, p);                 // Inicia a Matriz W de forma aleatória, como especificado//
    double** H = CriarMatriz(p, m);        // Cria a matriz H//
    double Er = 1, E_atual = 0, E_ant = 0;     // Define as variáveis dos erros, Erro relativo, Erro atual e Erro anterior, repectivamente//
    double** B   = CriarMatriz(n, m);      // Cria a Matriz B, que será a cópia de A//
    double** A_t = CriarMatriz(m, n);      // Cria a matriz A_t, que será a transposta de A//
    double** H_t = CriarMatriz(m, p);      // Cria a matriz H_t, que será a transposta de H//
    double** W_t = CriarMatriz(p, n);      // Cria a matriz W_t, que será a transposta de W//
    double** A_aprox = CriarMatriz(n, m);  // Cria a matriz A_aprox, que será a multiplicação de W*H//
    int iteracao = 0;                      // Define a variável iteracao, que conta quantas iterações foram feitas dentro do próximo loop//

    Copiar_Matriz(B, A, n, m);             // Faz o processo de copiar a matriz A em B//
    while((iteracao < 100) && (Er > 0.00001)){ // Loop principal da função, tem como parametros o iteracao e o erro relativo//
      E_atual = 0;      
     // printf("%d iteração\n\n", iteracao + 1);
      Inicializar_Zeros(H, p, m);          // Inicializando a matriz H//
      Inicializar_Zeros(W_t, p, n);        // Inicializando a Matriz transposta de W//
      // Estas inicializações acima são importantes para garantir que a cada iteração os valores da iteração anterior não influenciem no calculo da iteração atual //
      Copiar_Matriz(A, B, n, m);           // Restaura a matriz A para seus valores originais, a partir da cópia B //
      NormalizarColunas(W, n, p);          // Normaliza as colunas da matriz W//
      Decomposicao_QR_Multipla(W, n, p, m, H, A); // Realiza a decomposição QR multipla a partir de W e A, gerando a matriz resposta H //
      Maximo_0(H, p, m);                   // Torna a matriz H uma matriz não negativa//
      TransporMatriz(H, H_t, p, m);        // Torna a matriz H_t a transposta de H//
      TransporMatriz(B, A_t, n, m);        // Torna a matriz A_t a transposta de A//
      Decomposicao_QR_Multipla(H_t, m, p, n, W_t, A_t); // Novamente resolve a decomposição QR multipla de H_t e A_t para encontrar a matriz resposta W_t //
      TransporMatriz(W_t, W, p, n);        // Torna a matriz W em uma transposta da matriz W_t, encontrada anteriormente //
      Maximo_0(W, n, p);                   // Torna esta nova matriz W em uma matriz não negativa //
      Multiplicar_Matriz(W, H, A_aprox, n, m, p); // Faz a multiplicação de W por H, definindo a matriz A_aprox //
      for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++)
          E_atual += pow(B[i][j] - A_aprox[i][j], 2.0); // Define o erro atual como o quadrado da diferença entre a Matriz A original e a solução aproximada que chegamos nesta iteração //
      }
      Er = fabs((E_atual-E_ant));  // Erro relativo definido como a diferença entre o anterior e o atual, dividido pelo atual //
      E_ant = E_atual;                     // Erro anterior recebe o atual//
      //printf("Erro: %.5f\n\n", Er);
      iteracao += 1;                       // Prepara para a próxima interação//
    }
    DeletarMatriz(A_t, m);
    DeletarMatriz(H_t, m);                
    DeletarMatriz(W_t, p);                 
    DeletarMatriz(B, n);                  
    DeletarMatriz(A_aprox, n);            
    DeletarMatriz(H, p);                  
    return W;
}

void Segunda_tarefa(){
//----------------------------------------------------------------------------------------------------------------//
//------------------------------------------SEGUNDA TAREFA teste fornecido)---------------------------------------//
//----------------------------------------------------------------------------------------------------------------//
  double** A = CriarMatriz(3, 3);
  double** W = CriarMatriz(3, 2);

  A[0][0] = 0.3;
  A[1][0] = 0.5;
  A[2][0] = 0.4;
  A[0][1] = 0.6;
  A[1][1] = 0;
  A[2][1] = 0.8;
  A[0][2] = 0;
  A[1][2] = 1;
  A[2][2] = 0;

  // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e <= 36; e++){
    if (e == 18)
      printf("SEGUNDA TAREFA teste fornecido)");
    else
      printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  printf("\n");
  // fim do cabeçalho //

  W = Fatoracao_nao_Negativa(A, 3, 2, 3);
  printf("Matriz W: \n");
  printMatriz(W, 3, 2);

  printf("\n");

  DeletarMatriz(W, 3);
  DeletarMatriz(A, 3);

  printf("\n");
  printf("\n");

//----------------------------------------------------------------------------------------------------------------//
//-------------------------------------------SEGUNDA TAREFA teste criado)-----------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

  double** B = CriarMatriz(4, 4);
  double** M = CriarMatriz(4, 2);

  B[0][0] = 0.3;
  B[1][0] = 0.5;
  B[2][0] = 0.4;
  B[0][1] = 0.6;
  B[1][1] = 0;
  B[2][1] = 0.8;
  B[0][2] = 0;
  B[1][2] = 1;
  B[2][2] = 0;
  B[3][0] = 0.7;
  B[3][1] = 0;
  B[3][2] = 1;
  B[3][3] = 0.5;
  B[0][3] = 0.9;
  B[1][3] = 0;
  B[2][3] = 0.1;

  // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e <= 39; e++){
    if (e == 20)
      printf("SEGUNDA TAREFA teste criado)");
    else
      printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  printf("\n");
  // fim do cabeçalho //

  M = Fatoracao_nao_Negativa(B, 4, 2, 4);

  printf("Matriz M: \n");
  printMatriz(M, 4, 2);

  DeletarMatriz(B, 4);
  DeletarMatriz(M, 4);

  printf("\n");
  printf("\n");

}

//----------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------TAREFA FINAL---------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------//

// Esta é a função de leitura, responsável por ler os arquivos fornecidos com o enúnciado, contendo a base de dígitos//
// do MNIST, ela armazena os valores lidos na matriz A fornecida como parâmetro//
void Ler_arquivo(double** A, int n, int m, char *arquivo){
  FILE *fpointer = fopen(arquivo, "r");
  Inicializar_Zeros(A, n, m);
  char d = 'a';
  int q = 0;
  for(int k = 0; k < n; k++){
      for(int j = 0; d != '\r' && d != '\n' && d != EOF; j++ ){
        fscanf(fpointer,"%d%c", &q, &d);
        if(j < m)
          A[k][j] = (double)q/255.0;
      }
      d = '\0';
  }
  fclose(fpointer);
}

// A função abaixo é responsável por treinar os digitos com base nos aquivos fornecidos, ela realiza a fatoração //
// não negativa da matriz lida como a matriz de treino para o digito n, encontrando a matriz W correspondente e  //
// retornando esta. //

double** Treinar_matriz(int n, int m, int p, char *arquivo){
  double** A = CriarMatriz(n, m);
  double** W = CriarMatriz(n, p);

  Ler_arquivo(A, n, m, arquivo);
  W = Fatoracao_nao_Negativa(A, n, p, m);
  DeletarMatriz(A, n);

  return W;
}

// Por fim temos a função final deste exercício programa, a Reconhecer_digitos, que poderia ser também chamada de //
// Tarefa_principal ou correlacionados, pois é aqui que reconheceremos os digitos manúscritos e classificaremos a //
// eficiência do nosso código. //
void Reconhecer_digitos(int linhas, int n_treino, int n_teste, int amostras){
  int n = linhas;
  int ndig_treino = n_treino;
  int n_test = n_teste;
  int p = amostras;
  int i;
  double** W[10];      //Vetor de matrizes W, serão as matrizes treinadas. //
  for(i = 0; i < 10; i++){
    W[i] = CriarMatriz(ndig_treino, p);
    Inicializar_Zeros(W[i], ndig_treino, p);
  }
  double** A = CriarMatriz(n, n_test);
  double** B = CriarMatriz(n, n_test);

    // Cabeçalho //
  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");

  for (int e = 0; e <= 55; e++){
    if (e == 27)
      printf("TAREFA FINAL");
    else
      printf("-"); 
  }
  printf("\n");

  for (int e = 0; e < 67; e++)
    printf("-");
  printf("\n");
  printf("\n");
  // fim do cabeçalho //

  W[0] = Treinar_matriz(n, ndig_treino, p, "train_dig0.txt");
  W[1] = Treinar_matriz(n, ndig_treino, p, "train_dig1.txt");
  W[2] = Treinar_matriz(n, ndig_treino, p, "train_dig2.txt");
  W[3] = Treinar_matriz(n, ndig_treino, p, "train_dig3.txt");
  W[4] = Treinar_matriz(n, ndig_treino, p, "train_dig4.txt");
  W[5] = Treinar_matriz(n, ndig_treino, p, "train_dig5.txt");
  W[6] = Treinar_matriz(n, ndig_treino, p, "train_dig6.txt");
  W[7] = Treinar_matriz(n, ndig_treino, p, "train_dig7.txt");
  W[8] = Treinar_matriz(n, ndig_treino, p, "train_dig8.txt");
  W[9] = Treinar_matriz(n, ndig_treino, p, "train_dig9.txt");

  int *gabarito = (int*) malloc(n_test * sizeof(int));        // Vetor que armazena a resposta para o digito analisado//
  int *classificados = (int*) malloc(n_test * sizeof(int));   // Vetor que armazena o resultado para o digito analisado //
  double *erros = (double*) malloc(n_test * sizeof(double));  // Vetor que armazena os erros quadraticos //
  FILE*fpointer = fopen("test_index.txt", "r");               // Apontador de leitura //
  for(i = 0; i < n_test; i++){
    erros[i] = INFINITY;
    fscanf(fpointer, "%d", &gabarito[i]);
  }

  fclose(fpointer);

  double** H[10];                                   // Vetor para armazenar as 10 matrizes H que serão resolvidas//
  for(i = 0; i < 10; i++)
    H[i] = CriarMatriz(p, n_test);

  int k = 0, j = 0;
  Ler_arquivo(A, n, n_test, "test_images.txt");     // Processo de leitura da matriz a ser classificada //
  Copiar_Matriz(B, A, n, n_test);                   // Salvamos uma cópia desta matriz //
  double **WxH = CriarMatriz(n, n_test);            // Matriz que sera usada para encontrar o valor do erro //
  double Erro = 0, Erro2 = 0, somatoria = 0;

  for (i = 0; i < 10; i++){                         // Laço para testar cada digito possível //
    Copiar_Matriz(A, B, n, n_test);
    Decomposicao_QR_Multipla(W[i], n, p, n_test, H[i], A);
    Inicializar_Zeros(WxH, n, n_test); 

    Multiplicar_Matriz(W[i], H[i], WxH, n, n_test, p);

    for (k = 0; k < n_test; k++){
      somatoria = 0;
      for (j = 0; j < n; j++){
        Erro2 = A[j][k] - WxH[j][k];         // Erro ao quadrado é a diferença entre A e a matriz WxH que encontramos //
        somatoria += (Erro2*Erro2);
      }
      Erro = sqrt(somatoria);
      if (erros[k] > Erro){
        erros[k] = Erro;
        classificados[k] = i;
      }
    }
  }
  int acertos = 0;
  int digito_correto[10];
  int digito_total[10];

  for(i = 0; i < 10; i++){
    digito_correto[i] = 0;
    digito_total[i] = 0;
  }

  for(i = 0; i < n_test; i++){
    digito_total[gabarito[i]]++;
    if(gabarito[i] == classificados[i]){
      acertos++;
      digito_correto[gabarito[i]]++;
    }
  }
  
  
  for (i = 0; i < 10; i++)
    printf("           Digito %d: %6d/%5d corretos => %.3f%%\n", i, digito_correto[i], digito_total[i], (double)digito_correto[i]/digito_total[i] * 100);
  printf("              %6d/%5d corretos => %.3f%%\n\n", acertos, n_test, (double) acertos/n_test * 100);
  for(i = 0; i < 10; i++){
    DeletarMatriz(W[i], n);
    DeletarMatriz(H[i], p);
  }
  DeletarMatriz(A, n);
  DeletarMatriz(B, n);
  DeletarMatriz(WxH, n);

  free(gabarito);
  free(classificados);
  free(erros);

  printf("fim\n");
}

int main () {
  int n_test , ndig_treino, p;

  Primeira_Tarefa();
  Segunda_tarefa();

  printf("Digite o número de amostras por dígito: ");
  scanf("%d", &p);
  printf("Digite o número de digitos usados para treino: ");
  scanf("%d", &ndig_treino);
  printf("Digite o número de digitos a serem testados: ");
  scanf("%d", &n_test);

  printf("\n");

  Reconhecer_digitos(784, ndig_treino, n_test, p);

  return 0;

}