#coding:utf8
"""
Reconstruct part of the tree, may be for a node,

in decoding time, then derive k-best hypothesis list.

So in this theory, the parent of this part should be rule fetcher.


Input a list of tokens, then get a list of hyps,(targets,sites,costs)
sorted, pruned.

Original non-terminals lattice are using rules of MAINTOKEN -> hyps



- Raphael 2012.9
"""
import copy
import os,sys
from abraham.setting import setting
from abraham.logger import log_error,log

from gentile.model import GentileModel
from gentile.ruletable import GentileRuleTable
from gentile.model import GentileModel
from gentile.simplepruner import SimpleCubePruner
from chiropractor.hypothesis import Hypothesis
from gentile.ruletable import PSEUDO_COSTS

GLUECOST = PSEUDO_COSTS

class Reconstructor:
  """
  Reconstruct.
  """
  ruletable = None
  source = None
  areaTags = None
  model = None
  beamSize = None
  maxNTCount = None
  tag = None
  mapLattice = None
  mapSpanTag = None
  smode = None
  dependentAreas = None
  originalSites = None
  depraved = None

  area = None
  sourceString = None
  subTreeDistance = None


  def __init__(self, ruletable, model, sense, area, sourceString, subTreeDistance, hypStacks, source, areaTags, dependentAreas, depraved=False):
    """
    Init.
    """
    self.depraved = depraved
    self.source = source
    self.areaTags = areaTags
    # self.dependentAreas = dependentAreas
    self.originalSites = [idx for idx in range(len(source)) if isinstance(source[idx], tuple)]
    self.hypStacks = hypStacks
    self.ruletable = ruletable
    self.model = model
    # Reset smode.
    self.smode = model.smode
    model.smode = False
    self.beamSize = setting.size_beam
    self.maxNTCount = setting.reconstruction_max_nt
    self.mapLattice = {}
    self.mapSpanTag = {}
    self.sense = sense
    self.area = area
    self.sourceString = sourceString
    self.subTreeDistance = subTreeDistance
    self.tag = None

  def appendScoreToRules(self, rules):
    """
    get language cost and calculate score then append to rules.
    """
    newRules = []
    for rule in rules:
      translation, sites, costs = rule
      # Calculate lm score.
      words = [w for w in translation.split(" ") if not w.startswith("<<")]
      lmCost = self.model.getSentenseCost(words)
      # Calculate score.
      score = self.model.calculateScore(costs + [lmCost])
      newRules.append((translation, sites, costs, score))
    return newRules

  def sortStack(self, stack):
    """
    Sort a stack and prune with beam size, makes sure score are appended to
    each item of this stack.
    """
    stack.sort(key=lambda r:r[-2], reverse=True)
    return stack[:self.beamSize]

  def buildLexicalStack(self, token):
    """
    Build lexical (word) stack in the lattice.
    """
    sourceString = self.sense.tokens[token - 1][1]
    tag = self.sense.tokens[token - 1][0]
    rulesFound = self.ruletable.findBySource(sourceString, [], None)
    if not rulesFound:
      rulesFound = [self.ruletable.buildPsuedoRule(sourceString, [])]
    stack = self.appendScoreToRules(rulesFound)
    # Append empty site-to-nhyp dict.
    for nStack, rule in enumerate(stack):
      hyp = tuple(list(rule)+[{}])
      stack[nStack] = hyp
    stack = self.sortStack(stack)
    # Dirty way to allow DT words be translated to nothing.
    if tag == "DT":
      nullHyp = copy.deepcopy(stack[0])
      nullHyp = tuple([""]+list(nullHyp[1:]))
      stack.insert(0, nullHyp)
      stack = stack[:self.beamSize]

    return stack

  def findAllInferredCombinations(self, begin, width):
    """
    Calculate all possible inferred combinations.
    For example:
    Input (0,4), then get [(0,1,0),(1,1,0),(2,1,0),(3,1,0)],
    [(0,2,1)),(2,1,0),(3,1,0)]....

    The last 0 stands for treating as non-terminal,
    1 stands for treating as terminal.

    Mechanism:
    A example for 4 width, generate bits 0000 to 1110,

    in this list, 0 stands for terminals, a continious block of 1 stands for
    non-terminals.

    Just generate result for each items in this list.
    """
    combs = []
    for bits in range(0, 2**width-1):
      comb = []
      duringNT = False
      ntLength = 0
      ntBeginPos = 0
      ntCount = 0
      for iBit in range(0, width):
        if (bits >> iBit) % 2 == 0:
          # Check for Non-terminal end.
          if duringNT:
            comb.append((begin + ntBeginPos, ntLength, 1))
            ntLength = 0
            ntBeginPos = 0
            duringNT = False
            ntCount += 1
          # Append terminal.
          comb.append((begin + iBit, 1, 0))
        else:
          # Append Non-terminals.
          if not duringNT:
            ntBeginPos = iBit
            duringNT = True
          ntLength += 1
      # Check for Non-terminal end.
      if duringNT:
        comb.append((begin + ntBeginPos, ntLength, 1))
        ntCount += 1
      if ntCount <= self.maxNTCount:
        combs.append(comb)
    return combs

  def getSpanTag(self, span):
    """
    Calculate and cache tag for given span.

    @param span: (begin, width)
    """
    if span in self.mapSpanTag:
      return self.mapSpanTag[span]
    else:
      begin, width = span
      tokens = self.source[begin:begin+width]
      tokens = [t for t in tokens if isinstance(t, int)]
      tag = self.sense.getTokensTag(tokens)
      self.mapSpanTag[span] = tag
      return tag

  def getSpanWord(self, span):
    """
    Get words (actually word) for span.
    """
    tokenId = self.source[span[0]]
    assert tokenId > 0
    return self.sense.tokens[tokenId - 1][1]
    # if tokenId > 0:
    #   return self.sense.tokens[tokenId - 1][1]
    # else:
    #   linkedTokenId = self.sense.mapNodeToMainToken[-tokenId]
    #   return "[%s]" % self.sense.tokens[linkedTokenId - 1][0]

  def fetchRulesForComb(self, comb):
    """
    Convert input in format of [(begin, width, isNT), ...]
    to sourceString.
    Then try to fetch rules.
    If unable to find enough inferred support stacks,
    or no rule was found, just return (None,None)
    """
    patterns = []
    ntSpanList = []
    for spanEx in comb:
      begin, width, isNT = spanEx
      span = (begin, width)
      if isNT and span not in self.mapLattice:
        return None, None
      if isNT:
        tag = self.getSpanTag(span)
        patterns.append("[%s]" % tag)
        ntSpanList.append(span)
      else:
        # Original site could not be terminal
        if begin in self.originalSites:
          return None, None
        word = self.getSpanWord(span)
        patterns.append(word)
    sourceString = " ".join(patterns)
    rulesFound = self.ruletable.findBySource(sourceString, ntSpanList, None)
    if not rulesFound:
      return None, None
    # mapReplace = {}
    # sites = []
    # numX = 0
    # for idx, ntSpan in enumerate(ntSpanList):
    #   originalSite = self.originalSite(ntSpan)
    #   if originalSite:
    #     mapReplace["[X%d]" % idx] = "<<%d>>" % originalSite
    #   else:
    #     if idx != numX:
    #       mapReplace["[X%d]" % idx] = "[X%d]" % numX
    #     sites.append(ntSpan)
    #     numX += 1
    # # Repair all.
    # newRules = []
    # for rule in rulesFound:
    #   target, _, costs = rule
    #   for k in mapReplace:
    #     target = target.replace(k, mapReplace[k])
    #   newRules.append((target, sites, costs))
    return rulesFound, ntSpanList

  def buildSCFGStack(self, begin, width):
    """
    Build stack of non lexical lattice items.
    """
    stack = []
    # Find all combinations of all support or inferred items.
    combs = self.findAllInferredCombinations(begin, width)
    for comb in combs:
      rulesFound, ntSpanList = self.fetchRulesForComb(comb)
      if not rulesFound:
        continue
      # Merge products.
      simplePruner = SimpleCubePruner(self.model, rulesFound, ntSpanList, self.mapLattice)
      hyps = simplePruner.prune()
      # Prune with beam size.
      hyps = hyps[:self.beamSize]
      stack.extend(hyps)
    # If not found, then return None
    if not stack:
      return None
    # Only leave highest score of same translation.
    mapTranslationBestScore = {}
    for hyp in stack:
      translation = hyp[0]
      score = hyp[-2]
      if (translation not in mapTranslationBestScore or
          score > mapTranslationBestScore[translation]):
        mapTranslationBestScore[translation] = score
    idx = 0
    while idx < len(stack):
      hyp = stack[idx]
      if hyp[-2] < mapTranslationBestScore[hyp[0]]:
        stack.pop(idx)
      else:
        idx += 1
    # Sort and prune.
    return self.sortStack(stack)

  # def repairFinalRules(self, stack):
  #   """
  #   Restore a hyp stack to a list of rule, convert originalSite back to X.
  #   """
  #   sites = [-t for t in self.tokens if t < 0]
  #   newStack = []
  #   for hyp in stack:
  #     target, _, costs, _ = hyp
  #     for idx in range(len(sites)):
  #       target = target.replace("<<%d>>" % sites[idx], "[X%d]" % idx)
  #     newStack.append((target, sites, costs))
  #   return newStack

  def buildStackByGlueRule(self, supportLattice1, supportLattice2):
    """
    Build a pseudo lattice stack by given two support latice items.

    Output:
    [(
      "translation", [supportLattice1,supportLattice2],
      [COSTSUM - 10 * 6],
      score
    )]
    """
    # Check availibility of support lattice.
    if supportLattice1 not in self.mapLattice or supportLattice2 not in self.mapLattice:
      return None
    ntSpan = [supportLattice1, supportLattice2]
    translation = "[X0] [X1]"
    # Check if one of the two support lattice is original site.
    glueRule = (translation, ntSpan, GLUECOST)
    simplePruner = SimpleCubePruner(self.model, [glueRule], ntSpan, self.mapLattice)
    stack = simplePruner.prune()
    # Prune with beam size.
    stack = stack[:self.beamSize]
    return stack

  def putGlueRulesIntoWidthLayer(self, width):
    """
    Because no rules found in given width layer, put glue rules into it,
    and make CYK continue to run.
    """
    for begin in range(0, len(self.source) - width + 1):
      for lenLeftPart in range(1, width):
        stack = self.buildStackByGlueRule((begin, lenLeftPart), (begin + lenLeftPart, width - lenLeftPart))
        if stack:
          # print >> sys.stderr, (begin, lenLeftPart), (begin + lenLeftPart, width - lenLeftPart), len(stack)
          self.mapLattice[(begin, width)] = stack

  def convertHypsToLatticeStack(self, area):
    """
    Convert a list of hyps into lattice stack.
    """
    hyps = self.hypStacks[area]
    stack = []
    nHyp = 0
    for hyp in hyps:
      item = (hyp.translation, [], hyp.costs, hyp.score, {})
      stack.append(item)
      nHyp += 1
    return stack

  def convertLatticeStackToHyps(self, stack):
    """
    Convert a lattice stack to hyps.
    """
    hyps = []
    for item in stack:
      translation, _, costs, score, supportDict = item
      # TODO: save area sourceString, subTreeDistance
      hyp = Hypothesis(self.model, self.area, self.sourceString, None, None, self.subTreeDistance)
      hyp.sites = supportDict.keys()
      hyp.translation = translation
      hyp.translationTokens = translation.split(" ")
      hyp.costs = costs[:]
      # Build target and hypStack.
      target = translation
      # stackHypsSelected = {}
      # for linkedNode in supportDict:
      #   linkedHyp = self.hypStacks[linkedNode][supportDict[linkedNode]]
      #   stackHypsSelected[linkedNode] = linkedHyp
      #   subtranslation = linkedHyp.translation
      #   target = target.replace(subtranslation, "<<%d>>" % linkedNode)
      # hyp.target = target
      # hyp.stackHypsSelected = stackHypsSelected
      # Build score.
      hyp.buildScore()
      hyps.append(hyp)
    return hyps

  def depravedParsing(self):
    """
    Parse in depraved fashion.
    @return:
    """
    sites = []
    # Prepare base hypothesises
    for iPart in range(len(self.source)):
      part = self.source[part]
      if isinstance(part, tuple):
        self.mapLattice[(iPart, 1)] = self.convertHypsToLatticeStack(part)
      else:
        self.mapLattice[(iPart, 1)] = self.buildLexicalStack(part)
      sites.append((iPart, 1))


    psuedoRule = self.ruletable.buildPsuedoRule(None, sites)
    simplePruner = SimpleCubePruner(self.model, [psuedoRule], sites, self.mapLattice)
    simpleHyps = simplePruner.prune()
    hyps = self.convertLatticeStackToHyps(simpleHyps)
    return hyps


  def parse(self):
    """
    Parse converted token list, and finally derive a list of rules
    for given token list.
    tokens: [token id, ... ]
    context: tag
    """
    if self.depraved:
      return self.depravedParsing()
    for width in range(1, len(self.source) + 1):
      anyRulesFound = False
      if self.smode:
        # Maintain the smode.
        self.model.smode = (width == len(self.source))
      for idxBegin in range(0, len(self.source) - width + 1):
        if width == 1:
          token = self.source[idxBegin]
          if isinstance(token, tuple):
            # Non-terminal area
            # Put MainToken -> hyps into this lattice.
            assert token in self.hypStacks
            self.mapLattice[(idxBegin, 1)] = self.convertHypsToLatticeStack(token)
          else:
            # Single terminal token
            self.mapLattice[(idxBegin, 1)] = self.buildLexicalStack(token)
        else:
          stack = self.buildSCFGStack(idxBegin, width)
          if stack:
            anyRulesFound = True
            self.mapLattice[(idxBegin, width)] = stack
      # If no rules found in this width layer, inject glue rules.
      if width > 1 and not anyRulesFound:
        self.putGlueRulesIntoWidthLayer(width)
    if (0, len(self.source)) not in self.mapLattice:
      return []
    self.model.smode = self.smode
    hyps = self.convertLatticeStackToHyps(self.mapLattice[(0, len(self.source))])
    return hyps

def test():
  """
  Ekk.
  """
  ruletable = GentileRuleTable()
  model = GentileModel()
  from testdata.testtree import testtree
  def showRules(rules):
    if not rules:
      print rules
      return
    for r in rules:
      t,s,c = r
      print t,"|",s,c
    return ""
  print "--- buildLexicalStack ---"
  rc = Reconstructor(ruletable, model,testtree,[2,3,4,-7],"NN")
  print rc.tokens
  print rc.sense.tokens
  #print rc.buildLexicalStack(2)
  print "--- findAllInferredCombinations ---"
  #print rc.findAllInferredCombinations(1, 3)
  print "--- getSpanWord ---"
  #print rc.getSpanWord((1,1))
  print "--- getSpanTag ---"
  #print rc.getSpanTag((0,3))
  print "--- parse ---"
  print showRules(rc.parse())

if __name__ == "__main__":
  test()
