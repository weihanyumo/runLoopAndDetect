/*
* Filename: UpixBinaryTree.h
* Author:   Xiehong
* Date:     2012.9.29
*/

#ifndef _UPIX_BINARY_TREE_H_
#define _UPIX_BINARY_TREE_H_


#include "UpixAlgorithmTypeDef.h"


//#define INSERT_RECUSIVE



// Class of Binary Tree Node
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixBTNode
{
public:
	UpixBTNode(const T &_Data = T(), UpixBTNode *_pLeft = NULL, UpixBTNode *_pRight = NULL, UpixBTNode *_pParent = NULL)
	{
		Data = _Data;
		pParent = _pParent;
		pLeft = _pLeft;
		pRight = _pRight;
	}

	UpixBTNode(const UpixBTNode &other)
	{
		Data = other.Data;
		pParent = other.pParent;
		pLeft = other.pLeft;
		pRight = other.pRight;
	}

	UpixBTNode<T> &operator=(const UpixBTNode &rhs)
	{
		if (this != &rhs)
		{
			Data = rhs.Data;
			pParent = rhs.pParent;
			pLeft = rhs.pLeft;
			pRight = rhs.pRight;
		}
		return *this;
	}

	~UpixBTNode() {}

public:
	void SetLeftChild(UpixBTNode *pNode)
	{
		pLeft = pNode;
		if(pNode != NULL)
		{
			pNode->pParent = this;
		}
	}

	void SetRightChild(UpixBTNode *pNode)
	{
		pRight = pNode;
		if(pNode != NULL)
		{
			pNode->pParent = this;
		}
	}

	bool IsRoot() const { return pParent == NULL; }

	bool IsLeaf() const
	{
		return (pLeft == NULL) && (pRight == NULL);
	}

	bool IsLeftChild() const
	{
		return (pParent != NULL) && (pParent->pLeft == this);
	}

	bool IsRightChild() const
	{
		return (pParent != NULL) && (pParent->pRight == this);
	}

public:
	T Data;
	UpixBTNode *pParent;
	UpixBTNode *pLeft;
	UpixBTNode *pRight;
};




// Visitor
// true: traverse done and break
// false: continue
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixBTVisitor
{
public:
	UpixBTVisitor() {}
	virtual ~UpixBTVisitor() {}

	virtual bool Visit(UpixBTNode<T> *pNode) = 0;
};




// Class of Binary Tree
template <class T>
class UPIX_ALGORITHM_TEMPLATE UpixBinaryTree
{
public:
	UpixBinaryTree(UpixBTNode<T> *pRoot = NULL)
	{
		m_pRoot = pRoot;
	}

	~UpixBinaryTree()
	{
		Clear();
	}

public:
	void Clear()
	{
		if (m_pRoot != NULL)
		{
			Destory(m_pRoot);
			m_pRoot = NULL;
		}
	}

	UpixBTNode<T> *GetRoot() { return m_pRoot; }

	void Insert(const T &val)
	{
		UpixBTNode<T> *pNode = new UpixBTNode<T>(val);
		Insert(pNode);
	}

	void Insert(UpixBTNode<T> *pNode)
	{
		if (m_pRoot == NULL)
		{
			m_pRoot = pNode;
			return;
		}

#ifdef INSERT_RECUSIVE
		_Insert(m_pRoot, pNode);
#else
		_Insert(pNode);
#endif
	}

	bool Remove(const T &val)
	{
		UpixBTNode<T> *pNode = Find(val);
		return Remove(pNode);
	}

	bool Remove(UpixBTNode<T> *pNode)
	{
		if (pNode == NULL)
		{
			return false;
		}

		UpixBTNode<T> *pMaxNode = NULL;
		if (pNode->IsLeaf())
		{
			pMaxNode = NULL;
		}
		else if (pNode->pLeft != NULL && pNode->pRight == NULL)
		{
			pMaxNode = pNode->pLeft;
		}
		else if (pNode->pLeft == NULL && pNode->pRight != NULL)
		{
			pMaxNode = pNode->pRight;
		}
		else
		{
			pMaxNode = FindMax(pNode->pLeft);
			if(pMaxNode == pNode->pLeft)
			{
				pMaxNode->SetRightChild(pNode->pRight);
			}
			else
			{
				pMaxNode->pParent->SetRightChild(pMaxNode->pLeft);
				pMaxNode->SetLeftChild(pNode->pLeft);
				pMaxNode->SetRightChild(pNode->pRight);
			} 
		}

		// Link pMaxNode and pNode->pParent
		if(pNode->IsLeftChild())
		{
			pNode->pParent->SetLeftChild(pMaxNode);
		}
		else if (pNode->IsRightChild())
		{
			pNode->pParent->SetRightChild(pMaxNode);
		}
		else
		{
			m_pRoot = pMaxNode;
			if (m_pRoot != NULL)
			{
				m_pRoot->pParent = NULL;
			}
		}

		delete pNode;
		pNode = NULL;
		return true;
	}

	UpixBTNode<T> *Find(const T &val)
	{
		if (m_pRoot == NULL)
		{
			return NULL;
		}

		UpixBTNode<T> *pNode = m_pRoot;
		while (pNode != NULL)
		{
			if (val < pNode->Data)
			{
				pNode = pNode->pLeft;
			}
			else if (val > pNode->Data)
			{
				pNode = pNode->pRight;
			}
			else
			{
				return pNode;
			}
		}

		return NULL;
	}

	UpixBTNode<T> *FindMax(UpixBTNode<T> *pNode)
	{
		if (pNode == NULL)
		{
			return NULL;
		}

		UpixBTNode<T> *pMaxNode = pNode;
		while (pMaxNode->pRight != NULL)
		{
			pMaxNode = pMaxNode->pRight;
		}

		return pMaxNode;
	}

	UpixBTNode<T> *FindMin(UpixBTNode<T> *pNode)
	{
		if (pNode == NULL)
		{
			return NULL;
		}

		UpixBTNode<T> *pMinNode = pNode;
		while (pMinNode->pLeft != NULL)
		{
			pMinNode = pMinNode->pLeft;
		}

		return pMinNode;
	}

	int Depth(UpixBTNode<T> *pNode)
	{
		if (pNode != NULL)
		{
			int nDepthLeft = Depth(pNode->pLeft);
			int nDepthRight = Depth(pNode->pRight);
			int nDepth = (nDepthLeft > nDepthRight) ? nDepthLeft : nDepthRight;
			return nDepth + 1;
		}
		return 0;
	}

	int Depth() { return Depth(m_pRoot); }

	int LevelTraverse(UpixBTNode<T> *pNode)
	{
		int nLevel = 0;

		if (pNode == NULL)
		{
			return nLevel;
		}
		
		_Level(m_pRoot, pNode, nLevel, 1);
		return nLevel;
	}

	int Level(UpixBTNode<T> *pNode)
	{
		if (pNode == NULL)
		{
			return 0;
		}

		int nLevel = 1;
		UpixBTNode<T> *pParentNode = pNode;
		while (!pParentNode->IsRoot())
		{
			pParentNode = pParentNode->pParent;
			nLevel++;
		}

		if (pParentNode != m_pRoot)
		{
			nLevel = 0;
		}
		return nLevel;
	}

	int GetLeafCount(UpixBTNode<T> *pNode)
	{
		if (pNode != NULL)
		{
			if (pNode->IsLeaf())
			{
				return 1;
			}
			int nLeafLeft = GetLeafCount(pNode->pLeft);
			int nLeafRight = GetLeafCount(pNode->pRight);
			return nLeafLeft + nLeafRight;
		}
		return 0;
	}

	int GetLeafCount() { return GetLeafCount(m_pRoot); }

	// true: traverse done and break
	// false: continue
	bool PreOrderTraverse(UpixBTNode<T> *pNode, UpixBTVisitor<T> *pVisitor)
	{
		if (pNode != NULL)
		{
			if (pVisitor->Visit(pNode))
			{
				return true;
			}

			if (PreOrderTraverse(pNode->pLeft, pVisitor))
			{
				return true;
			}

			if (PreOrderTraverse(pNode->pRight, pVisitor))
			{
				return true;
			}
		}
		return false;
	}

	bool MidOrderTraverse(UpixBTNode<T> *pNode, UpixBTVisitor<T> *pVisitor)
	{
		if (pNode != NULL)
		{
			if (MidOrderTraverse(pNode->pLeft, pVisitor))
			{
				return true;
			}

			if (pVisitor->Visit(pNode))
			{
				return true;
			}

			if (MidOrderTraverse(pNode->pRight, pVisitor))
			{
				return true;
			}
		}
		return false;
	}

	bool PostOrderTraverse(UpixBTNode<T> *pNode, UpixBTVisitor<T> *pVisitor)
	{
		if (pNode != NULL)
		{
			if (PostOrderTraverse(pNode->pLeft, pVisitor))
			{
				return true;
			}

			if (PostOrderTraverse(pNode->pRight, pVisitor))
			{
				return true;
			}

			if (pVisitor->Visit(pNode))
			{
				return true;
			}
		}
		return false;
	}

	bool IsEmpty() { return m_pRoot == NULL; }

private:
	void Destory(UpixBTNode<T> *pNode)
	{
		if (pNode != NULL)
		{
			Destory(pNode->pLeft);
			Destory(pNode->pRight);
			delete pNode;
		}
	}

	void _Insert(UpixBTNode<T> *pNode, UpixBTNode<T> *pInsertNode)
	{
		if (pNode == NULL)
		{
			return;
		}

		if (pInsertNode->Data < pNode->Data)
		{
			if (pNode->pLeft == NULL)
			{
				pNode->SetLeftChild(pInsertNode);
				return;
			}
			else
			{
				_Insert(pNode->pLeft, pInsertNode);
			}
		}
		else if (pInsertNode->Data > pNode->Data)
		{
			if (pNode->pRight == NULL)
			{
				pNode->SetRightChild(pInsertNode);
				return;
			}
			else
			{
				_Insert(pNode->pRight, pInsertNode);
			}
		}
		else
		{
			return;
		}
	}

	void _Insert(UpixBTNode<T> *pNode)
	{
		UpixBTNode<T> *pCurNode = m_pRoot;
		while (pCurNode != NULL)
		{
			if (pNode->Data < pCurNode->Data)
			{
				if (pCurNode->pLeft != NULL)
				{
					pCurNode = pCurNode->pLeft;
				}
				else
				{
					pCurNode->SetLeftChild(pNode);
					return;
				}
			}
			else if (pNode->Data > pCurNode->Data)
			{
				if (pCurNode->pRight != NULL)
				{
					pCurNode = pCurNode->pRight;
				}
				else
				{
					pCurNode->SetRightChild(pNode);
					return;
				}
			}
			else
			{
				return;
			}
		}
	}

	bool _Level(UpixBTNode<T> *pNode, UpixBTNode<T> *pAimNode, int &nLevel, int nCurLevel)
	{
		if (pNode != NULL)
		{
			if (pNode->Data == pAimNode->Data)
			{
				nLevel = nCurLevel;
				return true;
			}

			nCurLevel++;

			if (_Level(pNode->pLeft, pAimNode, nLevel, nCurLevel))
			{
				return true;
			}

			if (_Level(pNode->pRight, pAimNode, nLevel, nCurLevel))
			{
				return true;
			}
		}

		return false;
	}

private:
	UpixBTNode<T> *m_pRoot;
};


#endif