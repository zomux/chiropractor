# coding: utf-8
"""
Fetch rules for a tree in Gentile.

For illustration, if a token list "o x o" exists.
"x" stands for non-terminal tokens, and "o" stands
for terminal tokens.

• Exactly Matching
  The source side of rule exactly matches the tokens list of the node,
  "o x o"

• Merged Matching
  Try to merge N levels of child nodes into the token list of this node,
  get "o o o", then use it as source side to search for rules.

• Reconstruct Matching(old method)
  If all above methods are failed, then try reconstruct matching,
  just build a new subtree for this node, to get chance of match a rule.
  Then this could be "o x o" --> "x x o" or "o x o" or "x x x".

• Reconstruct Matching(new method)
  CKY

• Depraved Matching
  If and only if all above methods are failed.
  Move all terminal tokens to separate child nodes and
  construct a pseudo rule "[tag-of-o] x [tag-of-o]" -> "x x x".

"""
import os,sys
from gentile.ruletable import GentileRuleTable
from chiropractor.hypothesis import Hypothesis
from abraham.setting import setting
from gentile.reconstructor import Reconstructor
import math
import itertools


class RuleFetcher:
  """
  Detemine pruning nodes of given tree
  Fetch all rules for each fragment derived by
  pruning nodes
  """

  ruletable = None
  model = None

  def __init__(self, ruletable, model):
    """
    Save the environment.

    @param sense: Sense Tree
    @type sense: SenseTree
    @param ruletable: ruletable object
    @type ruletable: GentileRuleTable
    """
    self.ruletable = ruletable
    self.model = model

  def findRulesBySourceString(self, sourceString, dependentAreas):
    """

    @param sourceString:
    @return:
    """
    rulesFound = self.ruletable.findBySource(sourceString, dependentAreas, None)
    return rulesFound
