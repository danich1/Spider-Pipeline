# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 16:38:16 2015

@author: Dave
"""

class KD_Tree:
    def __init__(self):
        self.root = None
        
    def insert(self,region,label):
        """
        Insert region into the KD-Tree

        region - a list of tuple containing start and end of a region
        label - the cassette or motif label
        """
        if self.root == None:
            self.root = Node(region,label)
        else:
           self.root.insert(region,label)
        leftsize = self.root.left.size() if self.root.left else 0
        rightsize = self.root.right.size() if self.root.right else 0
        if abs(rightsize-leftsize) > 1:
            self.resize()

    def resize(self):
        """
        Make the tree self balancing when adding new features
        """
        nodes = self.to_arr()
        node_index = int(len(nodes)/2)
        self.root = Node(nodes[node_index][0],nodes[node_index][1])
        nodes = nodes[:node_index] + nodes[node_index+1:]
        
        # reinsert nodes into KD tree
        for child in nodes:
            self.root.insert(child[0], child[1])

        
    def fits(self,query):
        """
        determine if a cassette can fit within the tree.
        """
        return self.root.fits(query) if self.root else True
        
    def get_label(self,index):
        """
        Return the node label

        index - query position
        """
        return self.root.get_label(index) if self.root else None
        
    def inorder(self):
        """
        Get the nodes in order
        """
        self.root.inorder([],True)
        
    def to_arr(self):
        """
        Display the KD Tree in array format
        """
        arr = []
        if self.root:
            self.root.inorder(arr)
        return arr
        
    def __str__(self):
        return str(self.root)
        
    def __repr__(self):
        return str(self.root)
        
    def size(self):
        """
        Get the size of the KD tree
        """
        return self.root.size()
        
class Node:
    def __init__(self,data,label):
        self.data = data
        self.left = None
        self.right = None
        self.label = label
        
    def insert(self,query,label):
        """
        Insert the child node into the KD tree.

        query - region to insert
        label - the name of the region to insert
        """
        if self.data[0] > query[1]:
            if self.left:
                self.left.insert(query,label)
            else:
                self.left = Node(query,label)
        elif self.data[1] < query[0]:
            if self.right:
                self.right.insert(query,label)
            else:
                self.right = Node(query,label)
        else:
            return None
            
    def get_label(self,index):
        """
        Get the name of the indicated node

        index - the region to index 
        """
        if self.data[0] <= index and self.data[1] >= index:
            return self.label
        elif self.data[0] > index:
            if self.left:
                return self.left.get_label(index)
        elif self.data[1] < index:
            if self.right:
                return self.right.get_label(index)
        return None
            
    def fits(self,query):
        """
        Determine if the node fits within the KD tree.

        query - the region to check if fits in tree
        """
        if self.data[0] > query[1]:
            if self.left:
                return self.left.fits(query)
            else:
                return True
        elif self.data[1] < query[0]:
            if self.right:
                return self.right.fits(query)
            else:
                return True
        else:
            return False
            
    def inorder(self,arr=[],verbose=False):
        """
        Get the nodes in order
    
        arr - the list to append
        verbse - flag to print out node labels
        """
        if self:
            if self.left:
                self.left.inorder(arr,verbose)
            if verbose:
                print(self)
            arr.append([list(self.data),self.label])
            if self.right:
                self.right.inorder(arr,verbose)

    def size(self):
        """
        Get the size of the tree.
        """
        leftSize = self.left.size() if self.left else 0
        rightSize = self.right.size() if self.right else 0
        return leftSize + 1 + rightSize
        
    def __str__(self):
        return str(self.data) + "," + self.label