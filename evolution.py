'''
# Bioinformatika - Boda Zsolt - EKKB9P
# Evlolúciós algoritmus az informatikában (gyakorlati kivitelezés)


'''

__author__ = 'Boda Zsolt'
__version__ = 1.0

from dataclasses import dataclass, field
from random import uniform, choice
from typing import List, Tuple, Union
from string import ascii_lowercase
from math import ceil
from uuid import uuid4, UUID
#from abc import ABC
#from functools import cached_property


# Ezekkel a változókkal lehet játszani és egész érdekes statisztikát kimutatni
# a geneárcióban részvevő példányok, a mutációs eslélye és hogy hány példányból
# generálunk új generációt arányából.
GENE_POOL_SIZE: int = 1000 # Ennyi példány lehet a generációban
GOAL_GENE_STRING: str = 'valami' # Ez a cél példányunk génjei
GENE_STRING_LEN: int = len(GOAL_GENE_STRING) # Gének száma (számoláshoz segít)
MUATION_RATE: float = 0.05 # Jelenleg 0.05 = 5%. Mutálodott gén esélye a gyerek képzésnél
# 0.3 = 30%. A generáció legjobb példányaiból, ennyi legjobbat választunk ki
GENE_POOL_PERCENTAGE: float = 0.3


@dataclass
class Instance:
    '''
        Példányt reprezentáló model. A példány tartalmaz egy egyedi azonosítót,
        amivel megkülönböztetjük a példányokat. Valamint géneket, amit egy string
        (karakter sorozat) reprezentál.
    '''
    genes: List[chr]
    _id: UUID = field(default_factory=uuid4)

    @staticmethod
    def generate_instance(gene_len: int):
        ''' Létrehoz egy teljesen új példányt, véletlenül választott génekkel. '''
        return Instance([choice(ascii_lowercase) for _ in range(gene_len)])

    def __eq__(self, o: object) -> bool:
        ''' Segítő metódus, megnézi hogy ugyanarról a kettő példányról van szó. '''
        if not isinstance(self, object):
            raise ValueError('Instance only can be compared to Instance!')

        return self._id == o._id

    def calc_fitness(self, instance) -> int:
        ''' Kiszámolja, hogy "milyen közel' áll az elérni kívánt példányhoz. '''
        if len(self.genes) != len(instance.genes): return -1

        return len([0 for (x, y) in zip(self.genes, instance.genes) if x == y])

    @staticmethod
    def create_child(parent_a, parent_b):
        ''' Létrehozza 2 "szülő" példány "gyermekét" / leszármazottát. '''
        genes: List[chr] = []

        for (a, b) in zip(parent_a.genes, parent_b.genes):
            chanche: float = uniform(0,1)

            if chanche <= MUATION_RATE:
                genes.append(choice(ascii_lowercase))
                continue

            gene_to_add: chr = a if chanche <= 0.5 else b
            genes.append(gene_to_add)

        return Instance(genes)



class GenePool:
    '''
        Generációt reprezentáló model. A generációkban példányok vannak.
        Ez a class moderálja az új generáció / genetikai lehetőségek létrehozását,
        valamint vizsgálja, hogy elértük már a "tökéletes géneket" / vég állapotot.
        Számontartja, hogy hány generáció leforgása alatt érjük el a kívánt eredméynt. 
    '''
    def __init__(
            self,
            goal: Instance,
            instance_factory: Instance,
            gene_len: int,
            mutation_rate: float,
            gene_pool_size: int,
            gene_pool_top: float
        ):
        self.goal = goal
        self.instance_factory = instance_factory
        self.gene_len = gene_len
        self.mutation_rate = mutation_rate
        self.gene_pool_size = gene_pool_size
        self.gene_pool_top = gene_pool_top

        self.cycle: int = 0
        self.reached_goal: bool = False
        self.current_instances: List[Instance] = self.init_instance_pool()

    def init_instance_pool(self) -> List[Instance]:
        ''' Elkészíti az első (nem különleges) generációt és annak példányait. '''
        return [self.instance_factory.generate_instance(self.gene_len) for _ in range(self.gene_pool_size)]

    def get_top_instances(self) -> List[Instance]:
        ''' Kiválasztja a legjobb példányokat a generációban, a megadott % mellett. '''
        number_of_top_instance: int = ceil(self.gene_pool_size * self.gene_pool_top)

        sorted_gene_pool: List[Instance] = sorted(
            self.current_instances,
            key=lambda instance: self.goal.calc_fitness(instance),
            reverse=True
            )
        top_instances: List[Instance] = sorted_gene_pool[:number_of_top_instance]
        return top_instances

    def get_parents(self, gene_pool: List[Instance]) -> Tuple[Instance, Instance]:
        ''' Generációból kiválaszt 2 szülő példányt, amik különbözőek. '''
        parent_a: Instance = choice(gene_pool)
        parent_b: Instance = choice(gene_pool)

        while parent_b == parent_a:
            parent_b = choice(gene_pool)

        return (parent_a, parent_b)

    def create_child(self, parent_a: Instance, parent_b: Instance) -> Instance:
        ''' Kiválasztott szülőkből készít gyermeket / leszármazottat. '''
        genes: List[chr] = []

        for (a, b) in zip(parent_a.genes, parent_b.genes):
            chanche: float = uniform(0,1)

            if chanche <= MUATION_RATE:
                genes.append(choice(ascii_lowercase))
                continue

            gene_to_add: chr = a if chanche <= 0.5 else b
            genes.append(gene_to_add)

        return Instance(genes)

    def get_best_instance(self, gene_pool: List[Instance]) -> Instance:
        ''' Kiválasztja a generáció "legjobb" példányát. '''
        sorted_gene_pool: List[Instance] = sorted(
            gene_pool, key=lambda instance: self.goal.calc_fitness(instance),
            reverse=True
            )
        return sorted_gene_pool[0]

    def run_cycle(self):
        ''' Generáció múlását szimulálja. '''
        self.cycle += 1
        print(f'Beginning cycle number #{self.cycle} ...')

        top_instances: List[Instance] = self.get_top_instances()

        new_instances: List[Instance] = []

        # Run until the new instances number equals to the gene pool size
        while len(new_instances) != self.gene_pool_size:
            # Get 2 different parents
            (parent_a, parent_b) = self.get_parents(top_instances)

            # create new instance
            child = self.create_child(parent_a, parent_b)

            # Check if new instance is goal
            if self.goal.calc_fitness(child) == GENE_STRING_LEN:
                print('\n\n', '-'*30, sep='')
                print(f'Reached goal in generation #{self.cycle},\nwith instance:\n{child}')
                print('\n\n')
                self.reached_goal = True
                break

            # Add child to new generation
            new_instances.append(child)
            
        if not self.reached_goal:
            # Display best instance
            best_instance: Instance = self.get_best_instance(new_instances)
            print(f'Best instance has the genes: {"".join(best_instance.genes)}')

        # Set current instances
        self.current_instances = new_instances


if __name__ == '__main__':
    goal_gene_string_to_list: List[chr] = [c for c in GOAL_GENE_STRING] 
    goal = Instance(goal_gene_string_to_list)
    pool_manager = GenePool(
        goal,
        Instance,
        GENE_STRING_LEN,
        MUATION_RATE,
        GENE_POOL_SIZE,
        GENE_POOL_PERCENTAGE
    )

    try:
        while not pool_manager.reached_goal:
            pool_manager.run_cycle()
    except KeyboardInterrupt as e:
        print('Stopping simulator...')