'''
Created on Oct. 9, 2023

@author: Matthew
'''

from . import SplineException
# from numbers import Number
from bisect import insort
from typing import Dict, Iterator, List, Optional, Self, Tuple

Number = int | float

class KnotVector():
    FIRST_NORMALIZED_OUTER_KNOT = 0
    LAST_NORMALIZED_OUTER_KNOT = 1
    
    def __init__(self,
                 inner_knots: Optional[List[Number] | Tuple[Number, ...] | Dict[Number, int] | Number] = None,
                 first_outer_knot_mul: int = 4,
                 last_outer_knot_mul: int = 4):
        self.first_outer_knot_mul = first_outer_knot_mul
        self.last_outer_knot_mul = last_outer_knot_mul
        if inner_knots is None: # None
            self.__unique_knots = []
            self.__multiplicities = {}
            self.__length = 0
        elif isinstance(inner_knots, Number): # Number
            self.__unique_knots = [inner_knots]
            self.__multiplicities = {inner_knots: 1}
            self.__length = 1
        elif isinstance(inner_knots, Dict): # Dict[Number, int]
            self.__unique_knots = list(inner_knots.keys())
            self.__unique_knots.sort()
            self.__multiplicities = inner_knots.copy()
            self.__length = sum(inner_knots.values())
        elif isinstance(inner_knots, List | Tuple): # List[Number] | Tuple[Number, ...]
            self.__unique_knots = list(set(inner_knots))
            self.__unique_knots.sort()
            self.__multiplicities = {unique_knot: 0 for unique_knot in self.__unique_knots}
            for unique_knot in inner_knots:
                self.__multiplicities[unique_knot] += 1
            self.__length = len(inner_knots)
        else:
            raise SplineException('KnotVector must contain at least one knot')
    
    @property
    def first_outer_knot_mul(self) -> int:
        return self.__first_outer_knot_mul
    
    @first_outer_knot_mul.setter
    def first_outer_knot_mul(self, value: int) -> None:
        if value < 1:
            raise SplineException('First outer knot multiplicity must be at least one')
        self.__first_outer_knot_mul = value
    
    @property
    def last_outer_knot_mul(self) -> int:
        return self.__last_outer_knot_mul
    
    @last_outer_knot_mul.setter
    def last_outer_knot_mul(self, value: int) -> None:
        if value < 1:
            raise SplineException('Last outer knot multiplicity must be at least one')
        self.__last_outer_knot_mul = value
    
    @property
    def inner_knot_total(self) -> int:
        '''
        Returns the total quantity of inner knots, excluding duplicates
        '''
        return len(self.__unique_knots)
    
    @property
    def total_unique_knots(self) -> int:
        '''
        Returns the total quantity of knots, inner and outer, excluding duplicates
        '''
        return self.inner_knot_total + 2
    
    @property
    def first_inner_knot(self) -> Number:
        '''
        The first unnormalized inner knot
        '''
        return self.__unique_knots[0]
    
    @property
    def last_inner_knot(self) -> Number:
        '''
        The last unnormalized inner knot
        '''
        return self.__unique_knots[-1]
    
    @property
    def first_unnormalized_outer_knot(self) -> Number:
        '''
        The first unnormalized outer knot
        '''
        return self.unnormalize(KnotVector.FIRST_NORMALIZED_OUTER_KNOT)
    
    @property
    def last_unnormalized_outer_knot(self) -> Number:
        '''
        The last unnormalized outer knot
        '''
        return self.unnormalize(KnotVector.LAST_NORMALIZED_OUTER_KNOT)
    
    def __len__(self):
        return self.first_outer_knot_mul + self.__length + self.last_outer_knot_mul
    
    def __iter__(self) -> Self:
        self.__iter_index = -1
        return self
    
    def __next__(self) -> Number:
        if self.__iter_index < 0:
            self.__iter_index = 0
            return KnotVector.FIRST_NORMALIZED_OUTER_KNOT
        elif self.__iter_index == self.inner_knot_total:
            self.__iter_index += 1
            return KnotVector.LAST_NORMALIZED_OUTER_KNOT
        elif self.__iter_index > self.inner_knot_total:
            raise StopIteration
        u = self.normalize(self.__unique_knots[self.__iter_index])
        self.__iter_index += 1
        return u
    
    def toUnnormalizedVector(self) -> List[Number]:
        '''
        Expands the given KnotVector into a list of unnormalized knots
        '''
        out = [self.first_unnormalized_outer_knot] * self.first_outer_knot_mul
        
        for unique_knot in self.__unique_knots:
            out.extend([unique_knot] * self.multiplicity(unique_knot))
        
        out.extend([self.last_unnormalized_outer_knot] * self.last_outer_knot_mul)
        return out
    
    def toNormalizedVector(self) -> List[Number]:
        '''
        Expands the given KnotVector into a list of normalized knots
        '''
        out = [KnotVector.FIRST_NORMALIZED_OUTER_KNOT] * self.first_outer_knot_mul
        
        for unique_knot in self.__unique_knots:
            out.extend([self.normalize(unique_knot)] * self.multiplicity(unique_knot))
        
        out.extend([KnotVector.LAST_NORMALIZED_OUTER_KNOT] * self.last_outer_knot_mul)
        return out
    
    def __getitem__(self, index: Number) -> Number:
        '''
        Returns the normalized knot of the specified index.
        Includes inner and outer knots.
        Negative indexing wraps around, similar to list indexing.
        Indexing out of bounds in either direction returns the nearest (outer) knot.
        '''
        if index < 0:
            index += len(self)
        elif index < self.first_outer_knot_mul:
            return KnotVector.FIRST_NORMALIZED_OUTER_KNOT
        
        index -= self.first_outer_knot_mul
        
        for unique_knot in self.__unique_knots:
            multiplicity = self.multiplicity(unique_knot)
            if 0 <= index < multiplicity:
                return self.normalize(unique_knot)
            index -= multiplicity
        
        return KnotVector.LAST_NORMALIZED_OUTER_KNOT
    
    def items(self) -> Iterator[Tuple[Number, int]]:
        return self.__multiplicities.items()
    
    def __contains__(self, unnormalized_knot: Number) -> bool:
        '''
        Checks if an unnormalized knot (inner or outer) exists in the knot vector
        '''
        return unnormalized_knot in self.__multiplicities \
            | {self.first_unnormalized_outer_knot, self.last_unnormalized_outer_knot}
    
    def getKnotIntervalIndices(self, u: Number) -> Optional[Tuple[int, int]]:
        '''
        Given some unnormalized knot u in the knot interval [u_k, u_k+1)
        and supposing that u_k is the i-th unique knot, returns (k, i).
        If u belongs to no knot interval, returns None
        '''
        L = self.first_outer_knot_mul + self.__length
        
        if u < self.first_unnormalized_outer_knot:
            return None
        elif u > self.last_inner_knot:
            return (len(self) - 1, self.total_unique_knots - 1)
        elif u < self.first_outer_knot_mul:
            return (self.first_outer_knot_mul - 1, 0)
        
        k, i = 0, 1
        for knot in self.__unique_knots:
            if knot > u:
                return k + self.first_outer_knot_mul - 1, i - 1
            elif knot == u:
                break
            elif k == L - 1:
                break
            k += self.multiplicity(knot)
            i += 1
        
        return k + self.first_outer_knot_mul, i
    
    def getUniqueKnots(self) -> List[Number]:
        '''Returns a copy of the internal list of unnormalized knots'''
        return [u for u in self.__unique_knots]
    
    def getMultiplicities(self) -> Dict[Number, int]:
        return {u: m for u, m in self.__multiplicities.items()}
    
    def multiplicity(self, knot: Number) -> int:
        '''Returns the multiplicity of the given unnormalized knot'''
        return self.__multiplicities[knot]
    
    # TODO Test insertknot
    def insertKnot(self,
                   knot: Number,
                   multiplicity: int = 1) \
                        -> int:
        '''
        Inserts or re-inserts the knot the specified quantity of times,
        thereby increasing its multiplicity
        
        Returns the knot's new multiplicity
        '''
        if multiplicity <= 0:
            raise SplineException('Inserting a knot requires specifying a multiplicity greater than 0')
        
        self.__length += multiplicity
        
        if knot in self.__multiplicities:
            self.__multiplicities[knot] += multiplicity
        else:
            self.__multiplicities[knot] = multiplicity
            insort(self.__unique_knots, knot)
        
        return self.multiplicity(knot)
    
    # TODO Test removeKnot
    def removeKnot(self,
                   knot: Number,
                   multiplicity: int = 1) \
                        -> int:
        '''
        Removes the knot the specified quantity of times,
        thereby reducing its multiplicity.
        Removes the knot entirely if its multiplicity is reduced to 0
        
        Returns the knot's new multiplicity
        '''
        if knot not in self.__multiplicities:
            raise SplineException("Can't remove a knot that doesn't exist (i.e, has multiplicity 0)")
        elif multiplicity <= 0:
            raise SplineException('Removing a knot requires specifying a multiplicity greater than 0')
        elif multiplicity > self.multiplicity(knot):
            raise SplineException("Can't reduce a knot's multiplicity to be less than 0")
        
        self.__length -= multiplicity
        
        if multiplicity == self.multiplicity(knot):
            self.__unique_knots.remove(knot)
            del self.__multiplicities[knot]
            return 0
        self.__multiplicities[knot] -= multiplicity
        return self.multiplicity(knot)
    
    def normalize(self, unnormalized_knot: Number) -> Number:
        '''
        Returns the normalized knot.
        In other words, the new knot u should satisfy
        KnotVector.FIRST_NORMALIZED_KNOT <= u <= KnotVector.LAST_NORMALIZED_KNOT
        '''
        return (unnormalized_knot - self.first_inner_knot + 1) / (self.last_inner_knot - self.first_inner_knot + 2)
    
    def unnormalize(self, normalized_knot: Number) -> Number:
        '''
        Reverses knot normalization, returning the unnormalized knot
        '''
        return normalized_knot * (self.last_inner_knot - self.first_inner_knot + 2) + self.first_inner_knot - 1
    
    def normalizeAll(self) -> None:
        '''
        Replaces the unnormalized knots with normalized knots.
        Any old normalized knots therefore become new unnormalized knots.
        Good for knot vectors that don't expect new knots to be inserted.
        '''
        u0, u1 = self.first_inner_knot, self.last_inner_knot
        denom = u1 - u0
        
        if denom == 0:
            raise ZeroDivisionError('Cannot normalize knot vector')
        
        self.__unique_knots = [(u - u0 + 1) / (u1 - u0 + 2) for u in self.__unique_knots]
        self.__multiplicities = {(u - u0 + 1) / (u1 - u0 + 2): m for u, m in self.__multiplicities.items()}
    
    def __str__(self) -> str:
        return str(self.toUnnormalizedVector())
    
    def __repr__(self) -> str:
        return f'KnotVector({self.__multiplicities})'
