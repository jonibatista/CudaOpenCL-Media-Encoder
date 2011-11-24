Cuda Media Encode
===================

A codificação de conteúdos multimédia (imagens, video e som) é determinante, pois permite reduzir o
tamanho dos conteúdos, o que reduz significativamente custos no armazenamento e transmissão dos
conteúdos. Contudo, a codificação de conteúdos acarreta custos computacionais significativos, tanto nas
codificação como na descodificação. Diga se que se tolera melhor um algoritmo de codificação que seja
Diga-se computacionalmente oneroso, mas cujo processo de descodificação seja relativemente leve. A razão deriva
do facto que um conteúdo ser codificado apenas uma vez (na criação), mas descodificado inúmeras vezes,
(por exemplo, um video do YouTube ou uma imagem JPEG/PNG/etc., codificados na criação e JPEG/PNG/etc.,
descodificados sempre que é solicitado o seu visionamento).
Neste projeto pretende-se adaptar uma implementação do algoritmo Vector Quantization (disponível
tar sobre a forma de código fonte C) C)(Vector Quantization) ao paradigma de programação manycore
computing, nomeadamente ao CUDA (etapa 1) e ao OpenCL (etapa 2).


Utilização
----------
	sh run *TIPO_DE_IMAGEM*

#### Tipos de imagem	
*all* - codifica todas as imagens

*small* - codifica a mais pequena

*medium* - codifica a de tamanho médio



Alterações Codificador
----------
	int **Image_orig;
	int **Image_out;
	int ysize[1], xsize[1];           /* The dimensions of the original image */

Eram globais e não havia necessidade

	void start_outputting_bits(); // removeu-se esta função desnecessária
