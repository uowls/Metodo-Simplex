Os códigos foram feitos como parte de um trabalho computacional para a matéria de MS428 da Unicamp para a graduação de matemática aplicada (2s 2025). 
O código chamado Simplex.c foi criado para uso em um exemplo específico e deveria ser desconsiderado na maioria dos casos, já que ele tem a pimeria fase do método simplex ignorada, e tem a partição basica inicial factivel hardcoded para um caso especifico (após a primeira fase estava já na solução otima, que não seria bom para o desenho qeu acabamos colocando no trabalho mostando como o método simplex caminha ao resolver).
O código principal se chama metodoSimplex.c, e contém o código que reslve o método simplex de duas fases, com uma implementação simples de um método para lidar com degeneração entre as fases.
metodoSimplex.h é o header e contém as definições das duas funções implementadas em metodoSimplex.c. 
Não foi implementada nenhuma função auxiliar, somente as funções de resolver e de receber problemas de PL pelo terminal (que também converte ele para a forma padrão)
O código exemplo.c é o que tem o main que chama o prompt resolvedor e foi usado para fazer os exemplos "normais".
