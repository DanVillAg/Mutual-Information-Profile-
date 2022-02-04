from sklearn.metrics import adjusted_mutual_info_score
import matplotlib.pyplot as plt
import random
import numpy as np
import seaborn as sns

#CHR1 con funciones programadas H(X) + H(Y) - H(X,Y)

def genera_conjuntos(seq, k, M=None):
    """
    genera los conjuntos para evaluar la información mutua
    seq: iterable con símbolos
    k: corrimiento
    M: tamaño de la ventana. Si es None entonces se toma una tercera parte de seq
    """
    if(M is None):
        M = int(len(seq)/3)
        
    C = []
    X, Y = [], []
    v = seq[:M]
    for x, y in zip(v, seq[k:k+M]):
        C.append((x,y))
        X.append(x)
        Y.append(y)
    return C, X, Y

def get_mutual_info(seq, K):
    '''
    Retorna solo el arrego de info mutua
    '''
    I = []
    #K = 300
    T = int(len(seq)/2)
    for k in range(K):
        C, X, Y = genera_conjuntos(seq, k, T)
        I.append(adjusted_mutual_info_score(X,Y))
    return I

def generate_ami_profile(seq, K = 300, name = 'Default'):
    '''
    grafica el arreglo con los valores para el perfil ami
    seq: iterable con símbolos
    k: corrimiento
    name: el nombre para la gráfica
    '''
    I = get_mutual_info(seq, K)

    plt.plot(I[5:])
    plt.savefig('/home/dan/Escritorio/GenomicaComputacional/Mutual-Information-Profile-/proyecto/figures/{0}.png'.format(name))
    plt.xlabel("k")
    plt.ylabel("Información Mutua")
    plt.tight_layout()
    # Mostrar la gráfica
    plt.show()

def multiple_ami(list_seq, K = 50, name = 'Default'):
    '''
    grafica el arreglo con los valores para el perfil ami
    list_seq: lista de secuencia, iterable con símbolos
    k: corrimiento
    name: el nombre para la gráfica
    '''

    for sequence in list_seq:
        I = []
        #K = 50
        print(len(sequence))
        T = int(len(sequence)/2)
        for k in range(K):
            C, X, Y = genera_conjuntos(sequence, k, T)
            I.append(adjusted_mutual_info_score(X,Y))

        plt.plot(I[5:])
    
    plt.title("Superposición de todos los cromosomas")
    plt.xlabel("k")
    plt.ylabel("Información Mutua")
    plt.tight_layout()
    plt.savefig('/home/dan/Escritorio/GenomicaComputacional/Mutual-Information-Profile-/proyecto/figures/{0}.png'.format(name))
    # Mostrar la gráfica
    plt.show()

def random_sampling_corr(seq, samples=1000000):
    '''
    retorna la correlación de x samples contra una mayoría de la cadena
    '''
    print(len(seq))
    subseq = seq
    L = int(len(subseq))
    mutual_info_samples = []
    print(L)
    print(subseq[0:500])
    for sample in range(0,100):
        interval = random.randrange(0, L-5000)
        small_sample = subseq[interval:interval+5001]
        I = []
        K = 50
        T = int(len(small_sample)/2)
        for k in range(K):
            C, X, Y = genera_conjuntos(seq, k, T)
            I.append(adjusted_mutual_info_score(X,Y))
        mutual_info_samples.append(I)
    mutual_info_samples = np.array(mutual_info_samples)
    print(mutual_info_samples.shape)
    print('ok :)')
    return mutual_info_samples

def get_corr_against_gen(sequence):
    '''
    Saca las correlaciones de cada muestra contra el total
    sequence:lista de secuencias
    '''
    names = ['genome_S_areus','genome_E_coli']
    for idx, seq in enumerate(sequence):
        seq_mutual_info = get_mutual_info(seq, K = 50)
        print('Doing random 100,000 sampling of 5k sequences of genome: {0}'.format(names[idx]))
        samples_mutual_info = random_sampling_corr(seq, samples=100000)

        corrs = []
        for sample in samples_mutual_info:
            corrs.append(np.correlate(seq_mutual_info,sample))
        
        corrs = np.array(corrs)
        sns.displot(corrs, kind="kde", rug=True)

        plt.title("Histograma de correlacion entre muestras y primer millon de bases")
        plt.xlabel("Coeficiente de Correlacion")
        plt.tight_layout()
        plt.savefig('/home/dan/Escritorio/GenomicaComputacional/Mutual-Information-Profile-/proyecto/figures/corrs_{0}.png'.format(names[idx]))
    # Mostrar la gráfica
    plt.show()

