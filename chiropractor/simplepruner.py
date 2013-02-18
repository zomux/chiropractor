"""
Cube Pruner , you knows

Simple means the rule is (target, sites, costs)
hypothesis is (translation, sites, costs, score)
"""
from abraham.setting import setting

from gentile.model import GentileModel
from chiropractor.hypothesis import tree_distance_probability

from abraham.logger import log_error,log



class SimpleCubePruner:
  """
  Class for running cube pruning for zion
  prune all hypothesises in one cube

  dimension description:
  - size should be count of sites + 1
  - got first dimension for rule dimension
  """
  model = None
  """@type: GentileModel"""
  # [ (fragment,sites,target,costs) , ... ]
  rules = None
  # [ [hyp,...], ...]
  stackHypothesis = None

  size = None
  # position [idx,...]
  # first dimension should be rules
  bestPosition = None
  # [(score,pasted position),...]
  heapEdges = None
  # [(pasted position,hyp),...]
  mapPositionHyp = None
  dimensionSize = 0
  needStop = False
  # expanded positions
  positionExpanded = None
  layerLimits = None
  layetLimitsLength = None
  
  def __init__(self,model,rules,sitesInvolved,stackHyps):
    """
    Initial the cube pruner
    @type model: GentileModel
    @type rules: list of (list,str,list)
    @type sitesInvolved: list of number
    @type stackHyps: dict of list of GentileHypothesis
    """
    self.model = model
    self.rules = rules
    self.stackHypothesis = stackHyps
    self.sitesInvolved = sitesInvolved
    self.size = setting.size_cube_pruning
    
    self.dimensionSize = len(self.sitesInvolved) + 1
    self.bestPosition=[0]*self.dimensionSize
    self.heapEdges = []
    self.mapPositionHyp = {}
    self.positionExpanded = []
    self.lastBestScore = None
    # self.layerLimits = linearlimits(len(rules), self.size)
    # self.layerLimitsLength = len(self.layerLimits)

  def getHypothesisByPosition(self, position):
    """
    Get Hypothesis by given position of cube
    1. First dimension stands for rules
    2. for each sites involved, if the site is the head node
    currently pruning, then we should find hyps in lexical stack.
    @type position: list of number
    @rtype: NovelDepStrHypothesis
    """
    stackBasicHyps = {}
    rule = self.rules[position[0]]
    sitesForRule = rule[1]
    for idxPos in range(1,len(position)):
      idxHyp = position[idxPos]
      site = self.sitesInvolved[idxPos-1]
      if site in sitesForRule:
          # Just get hyp in normal hyp list.
          stackBasicHyps[site] = self.stackHypothesis[site][idxHyp]

    hyp = self.mergeHypothesises(rule, stackBasicHyps)
    return hyp

  def mergeHypothesises(self, rule, stackHypsSelected):
    """
    Merge several hypothesises by rule and support hyps.
    And recalculate the score.
    """
    translation, sites, costs = rule
    costs = costs[:]
    gotEmptyTranslation = False
    branches, matchedBranches = 0, 0
    for isite in range(len(sites)):
      site = sites[isite]
      assert site in stackHypsSelected
      basehyp = stackHypsSelected[site]
      subtranslation = basehyp[0]
      if len(subtranslation) == 0:
        gotEmptyTranslation = True
      translation = translation.replace("[X%d]" % isite, subtranslation)
      # merge costs
      for icost in range(len(costs)):
        costs[icost] += basehyp[2][icost]
      # Merge tree distances
      subBranches, subMatchedBranches = basehyp[-1]
      branches += subBranches
      matchedBranches = subMatchedBranches
    treeDistanceCost = tree_distance_probability(branches, matchedBranches)
    costs[0] = treeDistanceCost
    words = [w for w in translation.split(" ")]
    lmCost = self.model.getSentenseCost(words)
    score = self.model.calculateScore(costs + [lmCost])

    return (translation, sites, costs, score, (branches, matchedBranches))

  def pastePosition(self,position):
    """
    return a fully number of a list of positions
    @type position: list of int
    @rtype: string
    """
    return ','.join('%d' % num for num in position)

  def unpastePosition(self,pastedpos):
    """
    return a fully number of a list of positions
    @type pastedpos: string
    @rtype: list of int
    """
    ret = [int(i) for i in pastedpos.split(",")]
    return ret
  
  def expandEdges(self):
    """
    expand edge lattices of cube
    save new hypothesises to self.heapEdges
    -------------
    Expand could be done by add each dimension by 1 , of course , if it can be
    """

    # first get all positions to expand
    # [pos,...]
    newedges = []
    for idxDimension in range(self.dimensionSize):
      curIndex = self.bestPosition[idxDimension]
      if idxDimension == 0:
        # rule dimension
        if curIndex+1 < len(self.rules):
          newedge = self.bestPosition[:]
          newedge[0] = curIndex+1
          newedges.append( newedge )
      else:
        # 1. check if greater than length of list
        # 2. check if this rule covers the site
        #    - if not covers, then expand edge in this dimension
        #      should have no means
        siteDimension = self.sitesInvolved[idxDimension-1]
        lenStack = len(self.stackHypothesis[siteDimension])
        siteCovered = siteDimension in self.rules[self.bestPosition[0]][1]
        if curIndex+1 < lenStack and siteCovered:
          newedge = self.bestPosition[:]
          newedge[idxDimension] = curIndex+1
          newedges.append( newedge )

    # for each new edge lattice , build hypothesis
    for position in newedges:
      pastedpos = self.pastePosition(position)
      if pastedpos in self.positionExpanded : continue
      #print position
      hyp = self.getHypothesisByPosition(position)
      # !!! to prevent one direction expanding
      # if hyp.getScore() == self.lastBestScore : continue
      self.mapPositionHyp[pastedpos] = hyp
      self.heapEdges.append( (hyp[-2],pastedpos) )
      self.positionExpanded.append(pastedpos)
    if len(self.heapEdges) == 0:
      # no new hyp could be found
      self.needStop = True

  def popBestHypothesis(self):
    """
    pop the best hypothesis from self.heapEdges and return it
    @rtype: NovelDepStrHypothesis
    """
    self.heapEdges.sort(reverse=True)
    score, pastedposOfBest = self.heapEdges.pop(0)
    self.bestPosition = self.unpastePosition(pastedposOfBest)
    # while True:
    #   score, pastedposOfBest = self.heapEdges.pop(0)
    #   self.bestPosition = self.unpastePosition(pastedposOfBest)
    #   iRule = self.bestPosition[0]
    #   if iRule < self.layerLimitsLength:
    #     if self.layerLimits[iRule] == 0 and len(self.heapEdges)>0:
    #       continue
    #     else:
    #       self.layerLimits[iRule] -= 1
    #       break
    #   else:
    #     break

    return self.mapPositionHyp[pastedposOfBest]

  def prune(self):
    """
    Function actually do the pruning process
    return a sorted hypothesis list of this fragment
    @rtype: list of NovelDepStrHypothesis
    """
    hypBest = self.getHypothesisByPosition(self.bestPosition)
    stackOutput = [hypBest]
    self.lastBestScore = hypBest[-2]

    while len(stackOutput) < self.size:
      # expand currently best lattice
      self.expandEdges()

      if self.needStop : break

      # pop the best hyp of edges to output
      hypBest = self.popBestHypothesis()
      self.lastBestScore = hypBest[-2]
      stackOutput.append( hypBest )

    stackOutput.sort( key=lambda x:x[-2] , reverse=True )
    return stackOutput