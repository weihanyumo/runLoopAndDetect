/*
* Filename: UpixMap.h
* Author:   Xiehong
* Date:     2013.6.22
*/

#ifndef _UPIX_MAP_H_
#define _UPIX_MAP_H_


#include "UpixAlgorithmTypeDef.h"



template <class KeyType, class ValType>
class UPIX_ALGORITHM_TEMPLATE UpixRBTNode
{
public:
	UpixRBTNode(const KeyType &key, const ValType &val) 
		: pLeftChild(NULL), pRightChild(NULL), pParent(NULL), Key(key), Value(val), bIsRedNode(true) 
	{}

	~UpixRBTNode() {}

public:
	void SetLeftChild(UpixRBTNode *pNode)
	{
		pLeftChild = pNode;
		if(pNode != NULL)
		{
			pNode->SetParent(this);
		}
	}

	void SetRightChild(UpixRBTNode *pNode)
	{
		pRightChild = pNode;
		if(pNode != NULL)
		{
			pNode->SetParent(this);
		}
	}

	void SetParent(UpixRBTNode *pNode) { pParent = pNode; }

	void SetValue(const ValType &val) { Value = val; }

	void SetRed() { bIsRedNode = true; }

	void SetBlack() { bIsRedNode = false; }

	UpixRBTNode *GetLeftChild() const { return pLeftChild; }

	UpixRBTNode *GetRightChild() const { return pRightChild; }

	UpixRBTNode *GetParent() const { return pParent; }

	ValType GetValue() const { return Value; }

	KeyType GetKey() const { return Key; }

	bool IsRoot() const { return pParent == NULL; }

	bool IsLeaf() const
	{
		return (pLeftChild == NULL) && (pRightChild == NULL);
	}

	bool IsLeftChild() const
	{
		return (pParent != NULL) && (pParent->GetLeftChild() == this);
	}

	bool IsRightChild() const
	{
		return (pParent != NULL) && (pParent->GetRightChild() == this);
	}

	int GetLevel() const
	{
		if(IsRoot())
		{
			return 1;
		}
		else
		{
			return GetParent()->GetLevel() + 1;
		}
	}

	bool IsRed() const
	{
		return bIsRedNode;
	}

	bool IsBlack() const
	{
		return !bIsRedNode;
	}

private:
	UpixRBTNode();

private:
	UpixRBTNode *pLeftChild;
	UpixRBTNode *pRightChild;
	UpixRBTNode *pParent;
	KeyType Key;
	ValType Value;
	bool bIsRedNode;
};



template <class KeyType, class ValType>
class UPIX_ALGORITHM_TEMPLATE UpixMap
{
public:
	UpixMap() : pRoot(NULL), nSize(0) {}

	~UpixMap() { Clear(); }

public:
	class Iterator
	{
	public:
		Iterator() : pRoot(NULL), pCur(NULL) {}

		Iterator(UpixRBTNode<KeyType, ValType> *root) : pRoot(root) { Reset(); }

		Iterator(const Iterator &src) : pRoot(src.pRoot), pCur(src.pCur) {}

		void Reset(bool atLowest = true)
		{
			if(atLowest)
			{
				pCur = GetMin(pRoot);
			}
			else
			{
				pCur = GetMax(pRoot);
			}
		}

		bool AtEnd() const
		{
			return pCur == NULL;
		}

		UpixRBTNode<KeyType, ValType> *GetNode()
		{
			return pCur;
		}

		Iterator &operator = (const Iterator &src)
		{
			pRoot = src.pRoot;
			pCur = src.pCur;
			return (*this);
		}

		void operator ++ (int) { Inc(); }

		void operator -- (int) { Dec(); }

		UpixRBTNode<KeyType, ValType> *operator -> () { return GetNode(); }

		UpixRBTNode<KeyType, ValType> &operator * ()
		{
			if(AtEnd())
			{
				throw "Iterator at end";
			}
			return *pCur;
		}

	private:
		UpixRBTNode<KeyType, ValType> *GetMin(UpixRBTNode<KeyType, ValType> *pNode)
		{
			while(pNode && pNode->GetLeftChild())
			{
				pNode = pNode->GetLeftChild();
			}
			return pNode;
		}

		UpixRBTNode<KeyType, ValType> *GetMax(UpixRBTNode<KeyType, ValType> *pNode)
		{
			while(pNode && pNode->GetRightChild())
			{
				pNode = pNode->GetRightChild();
			}
			return pNode;
		}

		void Inc()
		{
			if(pCur == NULL)
			{
				return;
			}

			if(pCur->GetRightChild() != NULL)
			{
				pCur = GetMin(pCur->GetRightChild());
			}
			else if(pCur->IsLeftChild())
			{
				pCur = pCur->GetParent();
			}
			else
			{
				while(pCur->IsRightChild())
				{
					pCur = pCur->GetParent();
				}
				pCur = pCur->GetParent();
			}
		}

		void Dec()
		{
			if(pCur == NULL)
			{
				return;
			}

			if(pCur->GetLeftChild() != NULL)
			{
				pCur = GetMax(pCur->GetLeftChild());
			}
			else if(pCur->IsRightChild())
			{
				pCur = pCur->GetParent();
			}
			else
			{
				while(pCur->IsLeftChild())
				{
					pCur = pCur->GetParent();
				}
				pCur = pCur->GetParent();
			}
		}

		UpixRBTNode<KeyType, ValType> *pRoot;
		UpixRBTNode<KeyType, ValType> *pCur;
	};

	class ParentFirstIterator
	{
	public:
		ParentFirstIterator() : pRoot(0), pCur(0) {}

		explicit ParentFirstIterator(UpixRBTNode<KeyType, ValType> *root) : pRoot(root), pCur(NULL) { Reset(); }

		void Reset() { pCur = pRoot; }

		bool AtEnd() const { return pCur == NULL; }

		UpixRBTNode<KeyType, ValType> *GetNode() { return pCur; }

		ParentFirstIterator &operator = (const ParentFirstIterator &src)
		{
			pRoot = src.pRoot;
			pCur = src.pCur;
			return (*this);
		}

		void operator ++ (int) { Inc(); }

		UpixRBTNode<KeyType, ValType> * operator -> () { return GetNode(); }

		UpixRBTNode<KeyType, ValType> & operator * ()
		{
			if(AtEnd())
			{
				throw "ParentFirstIterator at end";
			}
			return *GetNode();
		}

	private:
		void Inc()
		{
			if(pCur == NULL)
			{
				return;
			}

			if(pCur->GetLeftChild() != NULL)
			{
				pCur = pCur->GetLeftChild();
			}
			else if(pCur->GetRightChild() != NULL)
			{
				pCur = pCur->GetRightChild();
			}
			else
			{
				while(pCur != NULL)
				{
					if(pCur->IsLeftChild() && pCur->GetParent()->GetRightChild() != NULL)
					{
						pCur = pCur->GetParent()->GetRightChild();
						return;
					}
					pCur = pCur->GetParent();
				}
			}
		}

		UpixRBTNode<KeyType, ValType> *pRoot;
		UpixRBTNode<KeyType, ValType> *pCur;
	};
	
	class ParentLastIterator
	{
	public:
		ParentLastIterator() : pRoot(NULL), pCur(NULL) {}

		explicit ParentLastIterator(UpixRBTNode<KeyType, ValType> *root) : pRoot(root), pCur(NULL) { Reset(); }

		void Reset() { pCur = GetMin(pRoot); }

		bool AtEnd() const { return pCur == NULL; }

		UpixRBTNode<KeyType, ValType> *GetNode() { return pCur; }

		ParentLastIterator &operator = (const ParentLastIterator &src)
		{
			pRoot = src.pRoot;
			pCur = src.pCur;
			return (*this);
		}

		void operator ++ (int)
		{
			Inc();
		}

		UpixRBTNode<KeyType, ValType> * operator -> ()
		{
			return GetNode();
		}

		UpixRBTNode<KeyType, ValType> & operator * ()
		{
			if(AtEnd())
			{
				throw "ParentLastIterator at end";
			}
			return *GetNode();
		}

	private:
		UpixRBTNode<KeyType, ValType> *GetMin(UpixRBTNode<KeyType, ValType> *pNode)
		{
			while(pNode != NULL && (pNode->GetLeftChild() != NULL || pNode->GetRightChild() != NULL))
			{
				if(pNode->GetLeftChild() != NULL)
				{
					pNode = pNode->GetLeftChild();
				}
				else
				{
					pNode = pNode->GetRightChild();
				}
			}
			return pNode;
		}

		void Inc()
		{
			if(pCur == NULL)
			{
				return;
			}

			if(pCur->IsLeftChild() && pCur->GetParent()->GetRightChild() != NULL)
			{
				pCur = GetMin(pCur->GetParent()->GetRightChild());
			}
			else
			{
				pCur = pCur->GetParent();
			}
		}

		UpixRBTNode<KeyType, ValType> *pRoot;
		UpixRBTNode<KeyType, ValType> *pCur;
	};
	
	class AccessClass
	{
		friend class UpixMap<KeyType, ValType>;

	public:
		void operator = (const ValType &val)
		{
			Tree.Set(Key, val);
		}

		operator ValType()
		{
			UpixRBTNode<KeyType, ValType> *pNode = Tree.Find(Key);
			if(pNode == NULL)
			{
				throw "Item not found";
			}
			return pNode->GetValue();
		}

	private:
		AccessClass(UpixMap &tree, const KeyType &key) : Tree(tree), Key(key) {}
		AccessClass();

		UpixMap &Tree;
		const KeyType &Key;
	};

	bool Insert(const KeyType &keyNew, const ValType &val)
	{
		UpixRBTNode<KeyType, ValType> *pNewNode = new UpixRBTNode<KeyType, ValType>(keyNew, val);

		if(!Insert(pNewNode))
		{
			delete pNewNode;
			return false;
		}

		while(!pNewNode->IsRoot() && (pNewNode->GetParent()->IsRed()))
		{
			if(pNewNode->GetParent()->IsLeftChild())
			{
				UpixRBTNode<KeyType, ValType> *pNewNodesUncle = pNewNode->GetParent()->GetParent()->GetRightChild();
				if(pNewNodesUncle != 0 && pNewNodesUncle->IsRed())
				{
					pNewNode->GetParent()->SetBlack();
					pNewNodesUncle->SetBlack();
					pNewNode->GetParent()->GetParent()->SetRed();
					pNewNode = pNewNode->GetParent()->GetParent();
				}
				else
				{
					if(pNewNode->IsRightChild())
					{
						pNewNode = pNewNode->GetParent();
						RotateLeft(pNewNode);
					}
					pNewNode->GetParent()->SetBlack();
					pNewNode->GetParent()->GetParent()->SetRed();
					RotateRight(pNewNode->GetParent()->GetParent());
				}
			}
			else
			{
				UpixRBTNode<KeyType, ValType> *pNewNodesUncle = pNewNode->GetParent()->GetParent()->GetLeftChild();
				if(pNewNodesUncle != 0 && pNewNodesUncle->IsRed())
				{
					pNewNode->GetParent()->SetBlack();
					pNewNodesUncle->SetBlack();
					pNewNode->GetParent()->GetParent()->SetRed();
					pNewNode = pNewNode->GetParent()->GetParent();
				}
				else
				{
					if(pNewNode->IsLeftChild())
					{
						pNewNode = pNewNode->GetParent();
						RotateRight(pNewNode);
					}
					pNewNode->GetParent()->SetBlack();
					pNewNode->GetParent()->GetParent()->SetRed();
					RotateLeft(pNewNode->GetParent()->GetParent());
				}
			}
		}

		pRoot->SetBlack();
		return true;
	}

	void Set(const KeyType &key, const ValType &val)
	{
		UpixRBTNode<KeyType, ValType> *pNode = Find(key);
		if(pNode != NULL)
		{
			pNode->SetValue(val);
		}
		else
		{
			Insert(key, val);
		}
	}

	UpixRBTNode<KeyType, ValType> *Delink(const KeyType &key)
	{
		UpixRBTNode<KeyType, ValType> *pNode = Find(key);
		if(pNode == NULL)
		{
			return NULL;
		}

		while(pNode->GetRightChild())
		{
			RotateLeft(pNode);
		}

		UpixRBTNode<KeyType, ValType> *pLeft = pNode->GetLeftChild();
		if(pNode->IsLeftChild())
		{
			pNode->GetParent()->SetLeftChild(pLeft);
		}
		else if(pNode->IsRightChild())
		{
			pNode->GetParent()->SetRightChild(pLeft);
		}
		else
		{
			SetRoot(pLeft);
		}

		--nSize;
		return pNode;
	}

	bool Remove(const KeyType &key)
	{
		UpixRBTNode<KeyType, ValType> *pNode = Find(key);
		if(pNode == NULL)
		{
			return false;
		}

		while(pNode->GetRightChild() != NULL)
		{
			RotateLeft(pNode);
		}

		UpixRBTNode<KeyType, ValType> *pLeft = pNode->GetLeftChild();
		if(pNode->IsLeftChild())
		{
			pNode->GetParent()->SetLeftChild(pLeft);
		}
		else if(pNode->IsRightChild())
		{
			pNode->GetParent()->SetRightChild(pLeft);
		}
		else
		{
			SetRoot(pLeft);
		}

		delete pNode;
		--nSize;
		return true;
	}

	void Clear()
	{
		ParentLastIterator itr(GetParentLastIterator());
		while(!itr.AtEnd())
		{
			UpixRBTNode<KeyType, ValType> *pNode = itr.GetNode();
			itr++;
			delete pNode;
		}
		pRoot = NULL;
		nSize = 0;
	}

	bool IsEmpty() const { return pRoot == NULL; }

	UpixRBTNode<KeyType, ValType> *Find(const KeyType &keyToFind) const
	{
		UpixRBTNode<KeyType, ValType> *pNode = pRoot;
		while(pNode != NULL)
		{
			KeyType key(pNode->GetKey());

			if(keyToFind == key)
			{
				return pNode;
			}
			else if(keyToFind < key)
			{
				pNode = pNode->GetLeftChild();
			}
			else
			{
				pNode = pNode->GetRightChild();
			}
		}
		return NULL;
	}

	UpixRBTNode<KeyType, ValType> *GetRoot() const { return pRoot; }

	int GetSize() const { return nSize; }

	Iterator GetIterator()
	{
		Iterator itr(GetRoot());
		return itr;
	}

	ParentFirstIterator GetParentFirstIterator()
	{
		ParentFirstIterator itr(GetRoot());
		return itr;
	}

	ParentLastIterator GetParentLastIterator()
	{
		ParentLastIterator itr(GetRoot());
		return itr;
	}

	AccessClass operator [] (const KeyType &key)
	{
		return AccessClass(*this, key);
	}

private:
	explicit UpixMap(const UpixMap &src);

	UpixMap &operator = (const UpixMap &src);

	void SetRoot(UpixRBTNode<KeyType, ValType> *pNewRoot)
	{
		pRoot = pNewRoot;
		if(pRoot != NULL)
		{
			pRoot->SetParent(NULL);
			pRoot->SetBlack();
		}
	}

	bool Insert(UpixRBTNode<KeyType, ValType> *pNewNode)
	{
		bool result = true;
		if(pRoot == NULL)
		{
			SetRoot(pNewNode);
			nSize = 1;
		}
		else
		{
			UpixRBTNode<KeyType, ValType> *pNode = pRoot;
			KeyType keyNew = pNewNode->GetKey();
			while(pNode != NULL)
			{
				KeyType key(pNode->GetKey());
				if(keyNew == key)
				{
					result = false;
					pNode = NULL;
				}
				else if(keyNew < key)
				{
					if(pNode->GetLeftChild() == NULL)
					{
						pNode->SetLeftChild(pNewNode);
						pNode = NULL;
					}
					else
					{
						pNode = pNode->GetLeftChild();
					}
				}
				else
				{
					if(pNode->GetRightChild() == NULL)
					{
						pNode->SetRightChild(pNewNode);
						pNode = NULL;
					}
					else
					{
						pNode = pNode->GetRightChild();
					}
				}
			}

			if(result)
			{
				++nSize;
			}
		}
		return result;
	}

	void RotateLeft(UpixRBTNode<KeyType, ValType> *pNode)
	{
		UpixRBTNode<KeyType, ValType> *pRight = pNode->GetRightChild();
		pNode->SetRightChild(pRight->GetLeftChild());
		if(pNode->IsLeftChild())
		{
			pNode->GetParent()->SetLeftChild(pRight);
		}
		else if(pNode->IsRightChild())
		{
			pNode->GetParent()->SetRightChild(pRight);
		}
		else
		{
			SetRoot(pRight);
		}
		pRight->SetLeftChild(pNode);
	}

	void RotateRight(UpixRBTNode<KeyType, ValType> *pNode)
	{
		UpixRBTNode<KeyType, ValType> *pLeft = pNode->GetLeftChild();
		pNode->SetLeftChild(pLeft->GetRightChild());
		if(pNode->IsLeftChild())
		{
			pNode->GetParent()->SetLeftChild(pLeft);
		}
		else if(pNode->IsRightChild())
		{
			pNode->GetParent()->SetRightChild(pLeft);
		}
		else
		{
			SetRoot(pLeft);
		}
		pLeft->SetRightChild(pNode);
	}

private:
	UpixRBTNode<KeyType, ValType> *pRoot;
	int nSize;
};



#endif