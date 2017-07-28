
#ifndef ARRAY_H
#define ARRAY_H

#include <vector>

template<class T> class ConstArray	{

	public:
	virtual ~ConstArray() {}

	virtual int GetSize() const = 0;
	virtual const T& GetVal(int index) const = 0;
};

template<class T> class Array : public ConstArray<T>	{

	public:
	virtual ~Array() {}

	virtual T& operator[](int index) = 0;
};

template<class T> class HomogeneousArray : public ConstArray<T>	{

	public:
	HomogeneousArray(int insize, const T& invalue) : size(insize), value(invalue) {}
	~HomogeneousArray() {}

	int GetSize() const /*override*/ {return size;}
	const T& GetVal(int index) const /*override*/ {return value;}

	private:
	int size;
	const T& value;
};

template<class T> class SimpleArray : public Array<T>	{

	public:
	SimpleArray(int insize) : array(insize) {}
	virtual ~SimpleArray() {}

	int GetSize() const /*override*/ {return array.size();}
	T& operator[](int index) /*override*/ {return array[index];}
	const T& GetVal(int index) const /*override*/ {return array[index];}
	const vector<T>& GetArray() const {return array;}

	protected:
	vector<T> array;
};


template<class T> class ConstMixtureArray : public ConstArray<T>	{

	public:

	ConstMixtureArray(const Array<T>* incomponents, const Array<int>* inalloc) : components(incomponents), alloc(inalloc)	{
	}
	~ConstMixtureArray() {}

	int GetSize() const /*override*/ {return alloc->GetSize();}
	const T& GetVal(int i) const /*override*/ {
		return components->GetVal(alloc->GetVal(i));
	}

	private:
	const Array<T>* components;
	const Array<int>* alloc;
};

template<class T> class MixtureArray : public Array<T>	{

	public:

	MixtureArray(Array<T>* incomponents, const Array<int>* inalloc) : components(incomponents), alloc(inalloc)	{
	}
	~MixtureArray() {}

	int GetSize() const /*override*/ {return alloc->GetSize();}
	const T& GetVal(int i) const /*override*/ {
		return components->GetVal(alloc->GetVal(i));
	}
	T& operator[](int i) /*override*/ {return (*components)[alloc->GetVal(i)];}

	private:
	const Array<T>* components;
	const Array<int>* alloc;
};

#endif

