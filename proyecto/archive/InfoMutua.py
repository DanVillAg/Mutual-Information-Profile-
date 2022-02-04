import itertools
from collections import Counter, OrderedDict
import numpy as np
import re
from scipy.special import rel_entr, entr, kl_div


def k_spaced_bases(first, second, space, sequence):
  '''
  Dados dos nucleótidos y una separación entre ellos devuelve cuantos elementos
  así existen en la secuencia para poder calcular la entropía
  Params:
    :first: primer nucleótido para buscar
    :second: segundo nucleótido para buscar
    :space: la k-separacion entre first y second
    :sequence: la secuencia sobre la que se buscará
  Returns:
    el número de elementos de first-second que estan k-separados
  '''
  reg = r"" + first + ".{"+ str(space) +"}" + second
  result = re.findall(reg, sequence)
  return len(result)

def get_joint_probs(sequence,k):
  '''
  Dada un secuencia retorna un numpy.array bidimensional con las probs. 
  conjuntas p_k(X,Y) de las bases para las 16 combinaciones posibles:
      [
        [(AA) (AC) (AG) (AT)]
        [(CA) (CC) (CG) (CT)]
        [(GA) (GC) (GG) (GT)]
        [(TA) (TC) (TG) (TT)]
      ]
  Params:
    :sequence: la secuencia a revisar 
  Returns:
    :joint_probs: un arreglo 4x4 de las probabilidades conjuntas [[(AA)(AC)...(CA)...]]
  '''
  nucl = ['A', 'C', 'T', 'G']
  sum = 0
  nucl_comb = {}
  for subset in itertools.product(nucl,repeat=2):
      fi, sec = subset
      curr_num = k_spaced_bases(fi, sec, k, sequence)
      nucl_comb[subset] = curr_num
      sum += curr_num

  normalized_count = {k: v / sum for k, v in nucl_comb.items()}
  joint_probs = OrderedDict(sorted(normalized_count.items()))
  # print(joint_probs)
  joint_probs = np.array(list(joint_probs.values()))
  joint_probs = joint_probs/joint_probs.sum(axis=0,keepdims=1)
  joint_probs = np.reshape(joint_probs, (4,4))

  return joint_probs

def product_of_marginals(probs):
  '''
  Crea una matriz intermedia del producto de probs. marginales p(X)p(Y)
  para poder hacer el calculo de la info. mutua
  Params:
    :probs: vector de probabilidades de nucleotidos
  Returns:
    :prod_marginals: la matriz de probs. marginales p(X)p(Y) para nucleótidos
  '''
  alist = []
  nucl = ['A', 'C', 'T', 'G']
  for index, _ in enumerate(nucl):
    curr_nucl = probs[index]
    prod_probs = probs*curr_nucl
    alist.append(prod_probs)
  prod_marginals = np.array(alist)
  return prod_marginals

def sort_tuple_vals(my_tup):  
  my_tup.sort(key = lambda x: x[0])
  return my_tup

def generate_marginal_probs(sequence):
  '''
  Dada la secuencia genera la distribución de prob. de las 4 bases
  A C G T
  Params: 
    :sequence: str >> la secuencia a analizar
  Returns:
    :prob_dist: np.array >> la distribución marginal de las bases en la secuencia
  '''
  counts=Counter(sequence)

  prob_dist = sort_tuple_vals(counts.most_common())
  prob_dist= [prob for (nucl,prob) in prob_dist]
  prob_dist = np.array(prob_dist)
  prob_dist = prob_dist/len(seq)
  prob_dist = np.array(prob_dist)
  prob_dist = prob_dist/prob_dist.sum(axis=0,keepdims=1)
  return prob_dist


def mutual_info(sequence, k):
  '''
  Se obtiene la información mutua de las bases k-separadas en una 
  secuencia dada, usando la fórmula:

    I_k = Σ_x Σ_y p_k(X,Y)·log(p_k(X,Y)/p(X)p(Y))

  Params:
    :sequence: la secuencia
    :k: la k separación que se busca
  Returns:
    :mutual_info_k: la métrica de mut. info. para la k escogida
  '''
  marginals = generate_marginal_probs(sequence)
  marginal_product = product_of_marginals(marginals)
  joint_probs = get_joint_probs(sequence,k)

  log_base_4 = np.log(joint_probs/marginal_product) / np.log(4)
  print(joint_probs, 'joint\n')
  print(marginals, 'marginals \n')
  print(marginal_product, 'marg prod \n')
  print(log_base_4, 'div of log\n')
  print(np.multiply(joint_probs, log_base_4), 'mult\n')

  mutual_info_k = np.nansum(np.multiply(joint_probs, log_base_4))
  # rel_entr, entr, kl_div
  # mutual_info_k = entr(marginals) + entr(marginals) - rel_entr(marginals,marginals)
  return mutual_info_k

def generate_ami_profile(sequence, k=50):
  '''
  Dada una secuencia genera la info. mutua para un rango de k's
  Params: 
    :sequence: la secuencia a correr
    :k: la separación de las bases. 50  por default pero puede expandirse para ver el espectro
  Returns:
    :(k_steps,mu): tupla de numpy.array con valores para k y para info. mutua
  '''
  k_steps = np.array(list(range(0,k+1)))
  mu = []
  for k_val in k_steps:
    mu.append(mutual_info(sequence,k_val))
  mu = np.array(mu)

  return k_steps, mu
