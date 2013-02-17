import os
import sys

from abraham.treestruct import DepTreeStruct
from abraham.setting import setting
from gentile.tagconvertor import convert_tags_for_tokens

setting.load(["x_as_tag"])

PATTERN_SEPARATE_NTS = "NP NP,NP NP NP,NP VP,NP VBZ,NP NP VBZ,S VBZ,DT X".split(",")

class GeneralTree:
  """
  A class represents gernal tree structure.
  """
  nodes = None
  mapParent = None
  mapChildren = None
  root = None

  def __init__(self):
    """
    Initialize members.
    """
    self.nodes, self.mapParent, self.mapChildren = {}, {}, {}

  def node(self, id):
    """
    Get the node by given node id.
    """
    return self.nodes[id] if id in self.nodes else None

  def parent(self, id):
    """
    Get the id of parent node.
    """
    return self.mapParent[id] if id in self.mapParent else None

  def children(self, id):
    """
    Get the id of parent node.
    """
    return self.mapChildren[id] if id in self.mapChildren else []

  def buildChildrenMap(self):
    """
    Automatically build map of chilren by map of parents.
    """
    self.mapChildren = {}
    for nodeId in self.mapParent:
      parentId = self.mapParent[nodeId]
      self.mapChildren.setdefault(parentId, []).append(nodeId)

class PCFGTree:
  """
  A class represents PCFG tree.
  """
  tokens = None
  tree = None


  def __init__(self, textPCFG):
    """
    Initialize.
    """
    self.tokens = []
    self.tree = GeneralTree()
    self.parsePCFG(textPCFG)
    self.tree.buildChildrenMap()

  def parsePCFG(self, textPCFG):
    """
    Build the tree with PCFG result.
    """
    if textPCFG.startswith("( "):
      textPCFG = textPCFG[2:]
    if textPCFG.endswith(" )"):
      textPCFG = textPCFG[:-2]
    textPCFG = textPCFG.replace("(","( ").replace(")"," )").replace("  "," ")
    # Begin parsing.
    nodeId = 1
    depth = 0
    branchStack = []
    currentNode = None
    currentNodeId = None
    for part in textPCFG.split(" "):
      if part == "(":
        if not currentNode:
          continue
        currentNodeId = nodeId
        nodeId += 1
        self.tree.nodes[currentNodeId] = currentNode
        if branchStack:
          self.tree.mapParent[currentNodeId] = branchStack[-1]
        else:
          self.tree.root = currentNodeId
        branchStack.append(currentNodeId)
        # Clear stage.
        currentNodeId = None
        currentNode = None
      elif part == ")":
        if branchStack:
          currentNodeId = branchStack.pop()
      else:
        if not currentNode:
          currentNode = (part, None)
        else:
          tag, _ = currentNode
          # For terminal it should be (tag, token id)
          currentNode = (tag, len(self.tokens) + 1)
          self.tokens.append((tag, part))
          # Record Terminal Node.
          currentNodeId = nodeId
          nodeId += 1
          self.tree.nodes[currentNodeId] = currentNode
          if branchStack:
            self.tree.mapParent[currentNodeId] = branchStack[-1]
          else:
            self.tree.root = currentNodeId
          branchStack.append(currentNodeId)
          # Clear stage.
          currentNode = None

    

class SenseTree:
  """
  A class represents Sense tree.
  """
  tokens = None
  mapNodeToMainToken = None
  mapMainTokenToNode = None
  tree = None
  depTree = None
  xAsTag = None

  def __init__(self, textPCFG, textDep):
    """
    Initialize.
    """
    self.xAsTag = setting.x_as_tag
    self.parseSentence(textPCFG, textDep)

  def parseSentence(self, textPCFG, textDep):
    """
    Parse tree structure from a penn style PCFG tree and a dep tree.
    """
    pcfg = PCFGTree(textPCFG)
    self.tokens = pcfg.tokens
    # Check x as tag.
    if self.xAsTag:
      for idx, token in enumerate(self.tokens):
        self.tokens[idx] = ("X", token[1])

    cfgTree = pcfg.tree
    textTag = " ".join(["%s/%s" % t for t in self.tokens])
    depTree = DepTreeStruct(textTag, textDep)
    self.depTree = depTree
    leftUpTree = self.liftUpTerminals(cfgTree)

    # Test for liftUpTerminals
    # print "--- liftUpTerminals ---"
    # def tok(l):
    #   if type(l) != list:
    #     return l
    #   return " ".join([self.tokens[x-1][1] for x in l])
    # for nodeId in leftUpTree.mapChildren:
    #   print nodeId,tok(leftUpTree.node(nodeId)), ":", \
    #     map(tok, map(leftUpTree.node, leftUpTree.mapChildren[nodeId])),leftUpTree.mapChildren[nodeId]

    # assert len(depTree.nodes) == len(self.tokens) + 1
    self.tree = leftUpTree
    self.buildMainTokenMap(self.tree, depTree)
    # Test for buildMainTokenMap.
    # print " ".join(["%d:%s" % (k,self.tokens[self.mapNodeToMainToken[k]-1][1])
    #        for k in self.mapNodeToMainToken if self.mapNodeToMainToken[k]])

    self.killNonTerminals(self.tree, depTree)
    # print "--- killNonTerminals ---"
    # def tok(l):
    #   if type(l) != list:
    #     return l
    #   return " ".join([self.tokens[x-1][1] for x in l])
    # for nodeId in leftUpTree.mapChildren:
    #   print nodeId,tok(leftUpTree.node(nodeId)), ":", \
    #     map(tok, map(leftUpTree.node, leftUpTree.mapChildren[nodeId])),leftUpTree.mapChildren[nodeId]


  def buildMainTokenMap(self, tree, depTree, nodeId=None):
    """
    Build a map saves main token for each node by dependency relations.
    """
    if not nodeId:
      nodeId = tree.root
      self.mapNodeToMainToken = {}

    for childNodeId in tree.children(nodeId):
      self.buildMainTokenMap(tree, depTree, childNodeId)

    node = [t for t in tree.node(nodeId) if t > 0]
    if type(node) != list:
      # Just map it to none value.
      self.mapNodeToMainToken[nodeId] = None
    elif nodeId not in self.mapNodeToMainToken:
      # Count which token have most nodes depent on it.
      mapTokenDepentCount = {}
      for tokenId in node:
        mapTokenDepentCount.setdefault(tokenId, 0)
        dependentTokenId = depTree.getParentNode(tokenId)
        if dependentTokenId > 0 and dependentTokenId in node:
          mapTokenDepentCount.setdefault(dependentTokenId, 0)
          mapTokenDepentCount[dependentTokenId] += 1
      # Token wins with greatest count, and then greatest token id.
      mainTokenId = None
      maxCount = -1
      for tokenId in node:
        countOfToken = mapTokenDepentCount[tokenId]
        if countOfToken >= maxCount:
          mainTokenId = tokenId
          maxCount = countOfToken
      self.mapNodeToMainToken[nodeId] = mainTokenId

    if nodeId == tree.root:
      # Build mapMainTokenToNode.
      self.mapMainTokenToNode = {}
      for nodeId in self.mapNodeToMainToken:
        mainTokenId = self.mapNodeToMainToken[nodeId]
        if not mainTokenId:
          continue
        self.mapMainTokenToNode[mainTokenId] = nodeId

  def getTokensTag(self, tokens):
    """
    Get tag of tokens.

    @param tokens: [token id, ...]
    """
    if not tokens:
      return "X"
    mapTokenDepentCount = {}
    for tokenId in tokens:
      mapTokenDepentCount.setdefault(tokenId, 0)
      dependentTokenId = self.depTree.getParentNode(tokenId)
      if dependentTokenId > 0 and dependentTokenId in tokens:
        mapTokenDepentCount.setdefault(dependentTokenId, 0)
        mapTokenDepentCount[dependentTokenId] += 1
    # Token wins with greatest count, and then greatest token id.
    mainTokenId = None
    maxCount = -1
    for tokenId in tokens:
      countOfToken = mapTokenDepentCount[tokenId]
      if countOfToken >= maxCount:
        mainTokenId = tokenId
        maxCount = countOfToken
    return self.tokens[mainTokenId - 1][0]



  def killNonTerminals(self, tree, depTree, nodeId=None):
    """
    Remove all non-terminals in the tree by dependency relations.
    """
    if not nodeId:
      nodeId = tree.root
    childNodes = tree.children(nodeId)
    if not childNodes:
      return
    for childNodeId in childNodes:
      self.killNonTerminals(tree, depTree, childNodeId)

    node = tree.node(nodeId)
    if type(node) == list:
      return
    # The node with most depent count wins.
    mapNodeDepentCount = {}
    for childNodeId in childNodes:
      mapNodeDepentCount.setdefault(childNodeId, 0)
      mainTokenId = self.mapNodeToMainToken[childNodeId]
      depentTokenId = depTree.getParentNode(mainTokenId)
      if depentTokenId not in self.mapMainTokenToNode:
        continue
      depentNodeId = self.mapMainTokenToNode[depentTokenId]
      # Looply find the dependent node.
      while depentNodeId:
        if depentNodeId in childNodes:
          if depentNodeId != childNodeId:
            mapNodeDepentCount.setdefault(depentNodeId, 0)
            mapNodeDepentCount[depentNodeId] += 1
          break
        else:
          depentNodeId = tree.parent(depentNodeId)
    # Best node to replace non-terminal would be that with greatest depent count,
    # then the left.
    
    mainNodeId = None
    maxCount = -1
    for childNodeId in childNodes:
      countOfNode = mapNodeDepentCount[childNodeId]
      if countOfNode > maxCount:
        mainNodeId = childNodeId
        maxCount = countOfNode
    # if nodeId == 30:
    #   import pdb;pdb.set_trace()
    # Replace non-terminal node with that main node.
    # Set parent node of main node.
    if nodeId in tree.mapParent:
      parentId = tree.mapParent[nodeId]
      tree.mapParent[mainNodeId] = parentId
      childrenOfParent = tree.children(parentId)
      childrenOfParent[childrenOfParent.index(nodeId)] = mainNodeId
    else:
      # Current node is root.
      del tree.mapParent[mainNodeId]
    # Set parent node of child nodes.
    for childNodeId in childNodes:
      if childNodeId == mainNodeId:
        continue
      tree.mapParent[childNodeId] = mainNodeId
    # Set child nodes of main node.
    if len(childNodes) > 1:
      childNodes.remove(mainNodeId)
      childNodes.extend(tree.children(mainNodeId))
      childNodes.sort()
      tree.mapChildren[mainNodeId] = childNodes
    # Remove current non-terminal node.
    del tree.nodes[nodeId]
    if nodeId in tree.mapParent:
      del tree.mapParent[nodeId]
    if nodeId in tree.mapChildren:
      del tree.mapChildren[nodeId]
    # Modify tree root node.
    if nodeId == tree.root:
      tree.root = mainNodeId

  def mergeContinuousNTsForNode(self, nodeId):
    """
    Merge for given node.
    """
    tokens = self.tree.node(nodeId)
    continuousNTs = []
    currentNTs = []
    for token in tokens:
      if token > 0:
        if currentNTs:
          continuousNTs.append(currentNTs)
          currentNTs = []
      else:
        currentNTs.append(token)
    if currentNTs:
      continuousNTs.append(currentNTs)

    continuousNTs = [nts for nts in continuousNTs if len(nts)>1]

    for nts in continuousNTs:
      lastNodeId = -nts[-1]
      self.tree.nodes[lastNodeId] = nts[:-1] + self.tree.nodes[lastNodeId]
      reversedFrontNTs = nts[:-1]
      reversedFrontNTs.reverse()
      # Shrink tokens in node.
      self.tree.mapChildren.setdefault(lastNodeId, [])
      for nt in reversedFrontNTs:
        self.tree.nodes[nodeId].remove(nt)
        self.tree.mapChildren[lastNodeId].insert(0, -nt)

  def mergeContinuousNTs(self):
    """
    Merge continuous non-terminals into one non-terminal token,
    using deptree to decide the new tag,
    create a new node to cover these non-terminals.
    extract last linked node's tokens into this node. 
    """
    for nodeId in self.tree.nodes:
      self.mergeContinuousNTsForNode(nodeId)



  def upMergeNode(self, nodeId):
    """
    Up merge a node to its parent node.
    """
    if not self.tree.parent(nodeId):
      return
    # Change tokens in parent node.
    parentNode = self.tree.parent(nodeId)
    parentTokens = self.tree.node(parentNode)
    nodeTokens = self.tree.node(nodeId)
    idxLink = parentTokens.index(-nodeId)
    self.tree.nodes[parentNode] = parentTokens[:idxLink] + nodeTokens + parentTokens[idxLink+1:]
    # Change child nodes of parent node.
    self.tree.mapChildren[parentNode].remove(nodeId)
    self.tree.mapChildren[parentNode].extend(self.tree.children(nodeId))
    # Change parent node of child nodes.
    for childNodeId in self.tree.children(nodeId):
      self.tree.mapParent[childNodeId] = parentNode
    # Remove this node.
    del self.tree.nodes[nodeId]
    if nodeId in self.tree.mapParent:
      del self.tree.mapParent[nodeId]
    if nodeId in self.tree.mapChildren:
      del self.tree.mapChildren[nodeId]
    if nodeId in self.mapNodeToMainToken:
      mainToken = self.mapNodeToMainToken[nodeId]
      del self.mapNodeToMainToken[nodeId]
      del self.mapMainTokenToNode[mainToken]

  def upMergeAllConjNodes(self):
    """
    Merge all [IN] [TO] node.
    """
    nodesToMerge = []
    for nodeId in self.tree.nodes:
      mainToken = self.mapNodeToMainToken[nodeId]
      if self.tokens[mainToken-1][0] in ["TO", "IN"]:
        nodesToMerge.append(nodeId)
    for nodeId in nodesToMerge:
      self.upMergeNode(nodeId)

  def liftUpTerminals(self, tree):
    """
    Lift up terminals in a cfg tree to highest position they could reach.
    The node of tree to return should be a list of token id or just non-terminal symbol.
    So each should be '[1,2,3]' or just 'NP', and the token id starts at 1.
    """
    nodesToRemove = []
    for nodeId in tree.nodes:
      node = tree.nodes[nodeId]
      # Skip non-terminals.
      if type(node) == list or not node[1]:
        continue
      liftUpTo = tree.parent(nodeId)
      while True:
        # If the parent of this terminal node is an one branch node,
        # then lift it up more one level.
        if len(tree.children(liftUpTo)) == 1:
          parentOfLeftUp = tree.parent(liftUpTo)
          # Remove this node.
          nodesToRemove.append(liftUpTo)
          liftUpTo = parentOfLeftUp
        else:
          break
      if type(tree.nodes[liftUpTo]) != list:
        tree.nodes[liftUpTo] = []
      tree.nodes[liftUpTo].append(node[1])
      nodesToRemove.append(nodeId)
    for nodeId in nodesToRemove:
      if nodeId not in tree.nodes:
        continue
      del tree.nodes[nodeId]
      del tree.mapParent[nodeId]
    tree.buildChildrenMap()
    return tree

  def appendXToTree(self):
    """
    Append X(-nodeId) to corresponding position of each node.
    """
    tree = self.tree
    for nodeId in tree.nodes:
      childNodes = tree.children(nodeId)
      if not childNodes:
        continue
      # Insert -(child node id) to this node.
      node = tree.node(nodeId)
      assert type(node) == list
      for childNodeId in childNodes:
        mainTokenOfChild = self.mapNodeToMainToken[childNodeId]
        # Find position for insertion.
        pos = len(filter(lambda x: x<mainTokenOfChild, node))
        node.insert(pos, -childNodeId)

  def convertTags(self):
    """
    Convert tags for tokens.
    """
    self.tokens = convert_tags_for_tokens(self.tokens)

  def rebuildTopNode(self):
    """
    Rebuild top node, lift up the last period, and make sure the node not to be a pure non-terminal node.
    """
    rootId = self.tree.root
    rootTokens = self.tree.node(rootId)
    if self.tokens[rootTokens[-1] - 1][1] != ".":
      return
    periodId = rootTokens[-1]
    # Lift up the last period.
    newRootId = max(self.tree.nodes) + 1
    self.tree.nodes[newRootId] = [periodId]
    self.tree.mapParent[rootId] = newRootId
    self.tree.mapChildren[newRootId] = [rootId]
    self.mapNodeToMainToken[newRootId] = periodId
    self.tree.root = newRootId
    self.tree.nodes[rootId].remove(periodId)
    # Fix top node.
    if not self.tree.nodes[rootId]:
      self.tree.nodes[rootId] = "X"
      self.killNonTerminals(self.tree, self.depTree, rootId)

  def rebuildCommaNodes(self, nodeId=None):
    """
    Reconstruct nodes with comma.
    next , in NP , NP VBZ ---> next , IN [in NP] , VBZ [NP XXX]
    """
    if not nodeId:
      nodeId = self.tree.root

    for childNodeId in self.tree.children(nodeId):
      self.rebuildCommaNodes(childNodeId)

    tokens = self.tree.node(nodeId)
    commas = []
    for tokenId in tokens:
      if tokenId > 0 and self.tokens[tokenId - 1][1] == ",":
        commas.append(tokenId)
    if not commas:
      return
    # Get Partations.
    partations = []
    part = []
    for tokenId in tokens:
      if tokenId in commas:
        partations.append(part)
        part = []
      else:
        part.append(tokenId)
    if part:
      partations.append(part)
    # Create new nodes for parts with non-terminal.
    for idxPart in range(len(partations)):
      part = partations[idxPart]
      nts = [-t for t in part if t<0]
      if nts and len(part) > 1:
        # Create new node.
        newNodeId = max(self.tree.nodes) + 1
        self.tree.nodes[newNodeId] = part[:]
        self.tree.mapParent[newNodeId] = nodeId
        self.tree.mapChildren[newNodeId] = nts
        partations[idxPart] = [-newNodeId]
        if len(nts) == len(part):
          # The node is a pure non-terminal node.
          self.killNonTerminalsForXNode(newNodeId)
        self.buildMainTokenMap(self.tree, self.depTree, newNodeId)
        assert newNodeId in self.mapNodeToMainToken
    # Fix node.
    newTokens = []
    for i in range(len(partations)+len(commas)):
      if i % 2 == 0:
        newTokens.extend(partations[i/2])
      else:
        newTokens.extend([commas[(i-1)/2]])
    self.tree.nodes[nodeId] = newTokens
    self.tree.mapChildren[nodeId] = [-t for t in newTokens if t < 0]

  def killNonTerminalsForXNode(self, nodeId):
    """
    Remove all non-terminals in the tree by dependency relations.
    """
    tree = self.tree
    

    tokens = tree.node(nodeId)
    nts = [-t for t in tokens if t < 0]
    if len(nts) != len(tokens):
      return

    childNodes = tree.children(nodeId)
    # The node with most depent count wins.
    mapNodeDepentCount = {}
    for childNodeId in childNodes:
      mapNodeDepentCount.setdefault(childNodeId, 0)
      mainTokenId = self.mapNodeToMainToken[childNodeId]
      depentTokenId = self.depTree.getParentNode(mainTokenId)
      if depentTokenId not in self.mapMainTokenToNode:
        continue
      depentNodeId = self.mapMainTokenToNode[depentTokenId]
      # Looply find the dependent node.
      while depentNodeId:
        if depentNodeId in childNodes:
          if depentNodeId != childNodeId:
            mapNodeDepentCount.setdefault(depentNodeId, 0)
            mapNodeDepentCount[depentNodeId] += 1
          break
        else:
          depentNodeId = tree.parent(depentNodeId)
    # Best node to replace non-terminal would be that with greatest depent count,
    # then the left.
    
    mainNodeId = None
    maxCount = -1
    for childNodeId in childNodes:
      countOfNode = mapNodeDepentCount[childNodeId]
      if countOfNode > maxCount:
        mainNodeId = childNodeId
        maxCount = countOfNode
    # Replace non-terminal node with that main node.
    # Set parent node of main node.
    posMainNodeToken = tokens.index(-mainNodeId)
    newTokens = tokens[:posMainNodeToken]+tree.nodes[mainNodeId]+tokens[posMainNodeToken+1:]
    tree.nodes[nodeId] = newTokens
    tree.mapChildren[nodeId] = [-t for t in newTokens if t < 0]
    # Remove main node.
    del tree.nodes[mainNodeId]
    for childOfMainNode in tree.children(mainNodeId):
      tree.mapParent[childOfMainNode] = nodeId
    if mainNodeId in tree.mapChildren:
      del tree.mapChildren[mainNodeId]
    if mainNodeId in tree.mapParent:
      del tree.mapParent[mainNodeId]
    if mainNodeId in self.mapNodeToMainToken:
      del self.mapNodeToMainToken[mainNodeId]

  def separateContiniousNonTerminals(self, nodeId=None):
    """
    Separate continious non-terminals out from its node.
    Matches patterns:
    NP NP, NP NP NP, NP VP, NP VBZ,NP NP VBZ,S VBZ,DT X
    """
    if not nodeId:
      nodeId = self.tree.root

    for childNodeId in self.tree.children(nodeId):
      self.separateContiniousNonTerminals(childNodeId)

    tokens = self.tree.node(nodeId)
    nts = [-t for t in tokens if t < 0]
    if len(nts) < 2:
      return
    processed = False
    # Get Partations.
    partations = []
    part = []
    for tokenId in tokens:
      if tokenId > 0:
        if part:
          partations.append(part)
          part = []
      else:
        part.append(tokenId)
    if part:
      partations.append(part)
    # Create new nodes for parts with non-terminal.
    for idxPart in range(len(partations)):
      part = partations[idxPart]
      if len(part) <= 1:
        continue
      links = [-t for t in part]
      pattern = " ".join([self.tokens[self.mapNodeToMainToken[n]-1][0] for n in links])
      if pattern not in PATTERN_SEPARATE_NTS:
        continue
      # Create new node.
      newNodeId = max(self.tree.nodes) + 1
      self.tree.nodes[newNodeId] = part[:]
      self.tree.mapParent[newNodeId] = nodeId
      self.tree.mapChildren[newNodeId] = links
      
      # The node is a pure non-terminal node.
      self.killNonTerminalsForXNode(newNodeId)
      self.buildMainTokenMap(self.tree, self.depTree, newNodeId)
      assert newNodeId in self.mapNodeToMainToken
      # Fix node.
      tokens = self.tree.nodes[nodeId]
      posFirst = tokens.index(part[0])
      newTokens = tokens[:posFirst] + [-newNodeId] + tokens[posFirst+1:]
      for tokenId in part[1:]:
        newTokens.remove(tokenId)
      self.tree.nodes[nodeId] = newTokens
      processed = True
    if processed:
      self.tree.mapChildren[nodeId] = [-t for t in self.tree.nodes[nodeId] if t < 0]

  # def expandNodeByDep(self, nodeId):
  #   """
  #   Expand a node.
  #   """
    


  # def expandByDep(self):
  #   """
  #   Expand tree by dependencies.
  #   May be this could only be done for nodes with more than N terminals.
  #   """
  #   minExpandTerminals = setting.min_expand_terminals
  #   tree = self.tree
  #   for nodeId in tree.nodes:
  #     tokens = tree.node(nodeId)
  #     terminals = [t for t in tokens if t > 0]
  #     if len(terminals) > minExpandTerminals:
  #       self.expandNodeByDep(nodeId)

  def buildLevelMap(self, node=None, level=1):
    """
    Implementation for decoder.
    Build level map of sense tree.
    """
    if not node:
      self.maxLevel = 0
      self.mapLevelNodes = {}
      self.mapNodeLevel = {}
      node = self.tree.root
    if level > self.maxLevel:
      self.maxLevel = level
    self.mapLevelNodes.setdefault(level, []).append(node)
    self.mapNodeLevel[node] = level
    for childNode in self.tree.children(node):
      self.buildLevelMap(childNode, level + 1)

  def getMaxLevel(self):
    """
    Implementation for decoder.
    Return max level.
    """
    return self.maxLevel

  def getNodesByLevel(self, level):
    """
    Implementation for decoder.
    Return nodes of given level.
    """
    assert level <= self.maxLevel
    return self.mapLevelNodes[level]

  def getRootNode(self):
    """
    Get root node id.
    """
    return self.tree.root


def dump_sense_tree(senseTree):
  """
  Dump a sense tree to string.
  Format:
    nodeId | word-tokenId,word-tokenId,[node-mainTokenId]
    ...
  """
  # senseTree.appendXToTree()
  tree = senseTree.tree
  results = []
  mapTokenTag = {}
  for i,token in enumerate(senseTree.tokens):
    tag, word = token
    mapTokenTag[i+1] = tag
  for nodeId in tree.nodes:
    tokens = []
    for tokenId in tree.node(nodeId):
      if tokenId > 0:
        word = senseTree.tokens[tokenId-1][1]
        tokens.append("%s-%d" % (word, tokenId))
      else:
        linkedNodeId = -tokenId
        linkedTokenId = senseTree.mapNodeToMainToken[linkedNodeId]
        tokens.append("[%d-%d-%s]" % (linkedNodeId, linkedTokenId, mapTokenTag[linkedTokenId]))
    strTokens = " , ".join(tokens)
    results.append("%d | %s" % (nodeId, strTokens))
  return "\n".join(results)
  
  
if __name__ == '__main__':
  # Test mode.
  print "--- Test of PCFGTree ---"
  pcfg = "( (S (NP (JJ hexagonal) (NN hole) (NN 21j)) (VP (VBZ is) (VP (VBN formed) (PP (IN in) (NP (DT the) (NN end) (NN surface))) (PP (IN on) (NP (NP (DT the) (NN side)) (PP (IN of) (NP (JJ flange) (NN 21a))) (PP (IN in) (NP (NP (DT the) (NN center)) (PP (IN of) (NP (NP (PRP$ its) (, ,)) (PP (IN of) (NP (NN stud) (CD 21))))))))))) (. .)) )"
  dep = """amod(21j-3, hexagonal-1)
nn(21j-3, hole-2)
nsubjpass(formed-5, 21j-3)
auxpass(formed-5, is-4)
root(ROOT-0, formed-5)
prep(formed-5, in-6)
det(surface-9, the-7)
nn(surface-9, end-8)
pobj(in-6, surface-9)
prep(formed-5, on-10)
det(side-12, the-11)
pobj(on-10, side-12)
prep(side-12, of-13)
amod(21a-15, flange-14)
pobj(of-13, 21a-15)
prep(side-12, in-16)
det(center-18, the-17)
pobj(in-16, center-18)
prep(center-18, of-19)
poss(axis-21, its-20)
pobj(of-19, axis-21)
prep(axis-21, of-22)
pobj(of-22, stud-23)
num(stud-23, 21-24)"""
#   pcfg = "( (S (NP (JJ hexagonal) (NN hole) (NN 21j)) (VP (VBZ is) (VP (VBN formed) (PP (IN in) (DT the) ('' '') (NP (NN end)) ('' '')) (NP (NN surface)) (PP (IN on) (NP (NP (DT the) (NN side)) (PP (IN of) (NP (JJ flange) (NN 21a))) (PP (IN in) (NP (NP (DT the) (NN center)) (PP (IN of) (NP (NP (PRP$ its) (NN axis)) (PP (IN of) (NP (NN stud) (CD 21))))))))))) (. .)) )"
#   dep = """amod(21j-3, hexagonal-1)
# nn(21j-3, hole-2)
# nsubjpass(formed-5, 21j-3)
# auxpass(formed-5, is-4)
# root(ROOT-0, formed-5)
# prep(formed-5, in-6)
# dep(in-6, the-7)
# pobj(in-6, end-9)
# dobj(formed-5, surface-11)
# prep(formed-5, on-12)
# det(side-14, the-13)
# pobj(on-12, side-14)
# prep(side-14, of-15)
# amod(21a-17, flange-16)
# pobj(of-15, 21a-17)
# prep(side-14, in-18)
# det(center-20, the-19)
# pobj(in-18, center-20)
# prep(center-20, of-21)
# poss(axis-23, its-22)
# pobj(of-21, axis-23)
# prep(axis-23, of-24)
# pobj(of-24, stud-25)
# num(stud-25, 21-26)"""
  pcfg="( (S (NP (NP (DT the) (NN guide) (CD 5)) (, ,) (NP (NP (DT the) (NN guide) (NNS rollers) (CD 3)) (VP (VBN shown) (PP (IN in) (NP (NN fig))) (NP (NNP <DOT>) (CD 1))))) (VP (VBP are) (VP (VP (VBN designed) (S (VP (TO to) (VP (VB provide) (NP (JJ such) (NN guidance)))))) (PRN (-LRB- -LRB-) (CC and) (S (NP (NNS embodiments)) (ADVP (RB thereof)) (VP (MD will) (VP (VB be) (VP (VBN described) (ADVP (RB hereinunder)))))) (-RRB- -RRB-)))) (. .)) )"
  dep = """det(guide-2, the-1)
nsubj(designed-15, guide-2)
num(guide-2, 5-3)
cc(guide-2, and-4)
det(rollers-7, the-5)
nn(rollers-7, guide-6)
conj(guide-2, rollers-7)
num(rollers-7, 3-8)
partmod(rollers-7, shown-9)
prep(shown-9, in-10)
pobj(in-10, fig-11)
dobj(shown-9, <DOT>-12)
num(<DOT>-12, 1-13)
aux(designed-15, are-14)
root(ROOT-0, designed-15)
aux(provide-17, to-16)
xcomp(designed-15, provide-17)
amod(guidance-19, such-18)
dobj(provide-17, guidance-19)
dep(described-26, and-21)
nsubjpass(described-26, embodiments-22)
advmod(described-26, thereof-23)
aux(described-26, will-24)
auxpass(described-26, be-25)
parataxis(designed-15, described-26)
advmod(described-26, hereinunder-27)
"""
  # t = PCFGTree(pcfg)
  # print "tokens:"
  # print t.tokens
  # print "map:"
  # for nodeId in t.tree.mapChildren:
  #   print t.tree.node(nodeId), ":", \
  #     map(t.tree.node, t.tree.mapChildren[nodeId])
  # print "--- Test of SenseTree ---"
  t = SenseTree(pcfg,dep)
  t.rebuildTopNode()
  t.appendXToTree()
  t.upMergeAllConjNodes()
  t.rebuildCommaNodes()
  t.convertTags()
  t.separateContiniousNonTerminals()
  
  # t.mergeContinuousNTs()
  # print t.tree.nodes
  # print "--- getSpanTag ---"
  # print t.getTokensTag([1, 2, 3])
  # print "--- Test dump_sense_tree ---"
  print dump_sense_tree(t)