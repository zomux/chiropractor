# coding: utf8
"""
Class for present model and calculate the cost for rule and hypothesis,
This class knows the rule format, as Ruletable:
(target, sites, costs)

Implementation requires:

sortRules()
sortHypothesises()


"""
import os,sys

import math

from abraham.setting import setting
from abraham.logger import log_error,log
import platform
if platform.system() == "Darwin":
  from languagemodel.srilm import LanguageModel
else:
  from languagemodel.kenlmwrapper import LanguageModel


class Model:
  """
  class for presenting the model in Gentile Model
  """
  # { rule->cost object }
  # costs : static costs : dynamic costs : lm
  # score : costs * weights
  maxGram = None
  lm = None
  weights = None
  lenWeights = None
  weightsForRules = None
  cacheMode = False

  cacheLMProbs = None
  smode = False
  
  def __init__(self):
    """
    Initiate the language model
    initiate the model with weights in setting file
    """
    setting.load(["max_gram","file_lm","weights"])
    self.maxGram = int(setting.max_gram)
    sys.stderr.write("[GentileModel] Loading language model... \n")
    self.lm = LanguageModel(lm_file=setting.file_lm, n=setting.max_gram)
    self.weights = setting.weights
    self.lenWeights = len(self.weights)
    self.weightsForRules = self.weights[:-1]
    self.cacheLMProbs = {}

  def calculateRuleScore(self,rule):
    """
    calculate the cost of rule
    @type rule: string
    @rtype: object
    """
    costs = rule[2]
    return reduce(lambda x,y: x+y,[costs[i]*self.weightsForRules[i] for i in range(len(costs))])
    
  def sortRules(self,rules,limit=None):
    """
    For each [rules] , calculate the cost , and return sorted rules
    but only keep [limit] rules
    @type rules: list of string
    @type limit: number
    @rtype: list of string
    """
    if not limit:
      limit = setting.size_beam
      
    list_scores = []
    list_indexes = range(len(rules))

    for rule in rules:
      list_scores.append(self.calculateRuleScore(rule))


    list_indexes.sort(key=lambda x:list_scores[x],reverse=True)
    return [rules[i] for i in list_indexes[:limit]]

  def sortHypothesises(self, hyps, limit=None):
    """
    sort hypothesises

    @type hyps: list of GentileHypothesises
    """
    if not limit:
      limit = setting.size_cube_pruning
      
    hyps.sort(key=lambda x: x.score,reverse=True)
    return hyps[:limit]

  def calculateHypothesisScore(self,hyp):
    """
    Get score of a list of costs by plus together with weights

    @type hyp: GentileHypothesis
    """
    return reduce(lambda x,y: x+y,[hyp.costs[i]*self.weights[i] for i in range(self.lenWeights)])
  
  def calculateScore(self, costs):
    """
    Calculate score of a pure cost list.
    """
    return sum([costs[i]*self.weights[i] for i in range(len(costs))])

  def getSentenseCost(self,tokens):
    """
    @type tokens: list of string
    @rtype: float
    """
    if not tokens:
      return -0.15
    
    if self.smode:
      tokens.insert(0, "<s>")
      tokens.append("</s>")

    if len(tokens) == 1 and tokens[0] == "":
      return -0.15

    return self.lm.tokensProbability(tokens)

    # if len(tokens) < self.maxGram:
    #   iters = 1
    # else:
    #   iters = len(tokens) - self.maxGram + 1
    
    # prob = 0.0

    # for ibegin in range(iters):
    #   words = tokens[ibegin:ibegin+self.maxGram]
    #   if self.cacheMode:
    #     hashWords = hash(" ".join(words))
    #     try:
    #       prob_words = self.cacheLMProbs[hashWords]
    #     except :
    #       # could not find in cache
    #       prob_words = self.lm.readNGram(words)
    #       self.cacheLMProbs[hashWords] = prob_words
    #     prob += prob_words
    #   else:
    #     prob += self.lm.readNGram(words)
    # return prob

  # def getSentenseAverageCost(self,tokens):
  #   """
  #   @type tokens: list of string
  #   @rtype: float
  #   """
  #   if len(tokens) == 1 and tokens[0] == "":
  #     return -0.15
      
  #   if len(tokens) < self.maxGram:
  #     iters = 1
  #   else:
  #     iters = len(tokens) - self.maxGram + 1
    
  #   prob = 0.0

  #   for ibegin in range(iters):
  #     words = tokens[ibegin:ibegin+self.maxGram]
  #     if self.cacheMode:
  #       hashWords = hash(" ".join(words))
  #       try:
  #         prob_words = self.cacheLMProbs[hashWords]
  #       except :
  #         # could not find in cache
  #         prob_words = self.lm.readNGram(words)
  #         self.cacheLMProbs[hashWords] = prob_words
  #       prob += prob_words
  #     else:
  #       prob += self.lm.readNGram(words)
  #   return prob/iters

  def clearCache(self):
    """
    clear lm cache
    !!! deprecated
    """
    pass
    # self.cacheLMProbs.clear()