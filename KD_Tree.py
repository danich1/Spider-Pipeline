# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 16:38:16 2015

@author: Dave
"""
class KD_Tree:
    def __init__(self):
        self.root = None
        
    def insert(self,region,label):
        if self.root == None:
            self.root = Node(region,label)
        else:
           self.root.insert(region,label)
        leftsize = self.root.left.size() if self.root.left else 0
        rightsize = self.root.right.size() if self.root.right else 0
        if abs(rightsize-leftsize) > 1:
            self.resize()

    def resize(self):
        nodes = self.to_arr()
        self.root = Node(nodes[len(nodes)/2][0],nodes[len(nodes)/2][1])
        nodes = nodes[:len(nodes)/2] + nodes[len(nodes)/2+1:]
        map(self.root.insert,map(lambda x: x[0],nodes),map(lambda x: x[1],nodes))
        
    def fits(self,query):
        return self.root.fits(query) if self.root else True
        
    def get_label(self,index):
        return self.root.get_label(index) if self.root else None
        
    def inorder(self):
        self.root.inorder([],True)
        
    def to_arr(self):
        arr = []
        if self.root:
            self.root.inorder(arr)
        return arr
        
    def __str__(self):
        return str(self.root)
        
    def __repr__(self):
        return str(self.root)
        
    def size(self):
        return self.root.size()
        
class Node:
    def __init__(self,data,label):
        self.data = data
        self.left = None
        self.right = None
        self.label = label
        
    def insert(self,query,label):
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
        if self:
            if self.left:
                self.left.inorder(arr,verbose)
            if verbose:
                print(self)
            arr.append([list(self.data),self.label])
            if self.right:
                self.right.inorder(arr,verbose)

    def size(self):
        leftSize = self.left.size() if self.left else 0
        rightSize = self.right.size() if self.right else 0
        return leftSize + 1 + rightSize
        
    def __str__(self):
        return str(self.data) + "," + self.label