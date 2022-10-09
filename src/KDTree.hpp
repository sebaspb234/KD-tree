// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"
using namespace std;

template <size_t N, typename ElemType>
struct KDTreeNode
{
    Point<N> point; //Punto del nodo de N dimensiones
    ElemType value;
    KDTreeNode<N, ElemType>* nodes[2];
    KDTreeNode(Point<N> p, ElemType x)
    {
        point = p;
        value = x;
        nodes[0] = nodes[1] = nullptr;
    }
};


template <size_t N, typename ElemType>
class KDTree {
public:
    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();

    ~KDTree();

    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);

    size_t dimension() const;

    size_t size() const;
    bool empty() const;

    bool contains(const Point<N>& pt) const;

    void insert(const Point<N>& pt, const ElemType& value);

    ElemType& operator[](const Point<N>& pt);

    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;

    ElemType knn_value(const Point<N>& key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N>& key, size_t k) const;

    bool find(const Point<N>& pt, KDTreeNode<N, ElemType>**& p);

    mutable KDTreeNode<N, ElemType>* root;
    void neighbors(const Point<N> key, KDTreeNode<N, ElemType>* current_node, vector<pair<double, KDTreeNode<N, ElemType>*>>& nearest_neighbors_candidates, int depth, int k) const;
    KDTreeNode<N, ElemType>* copy(KDTreeNode<N, ElemType>* r);
    void deleteTree(KDTreeNode<N, ElemType>* r);
private:
    size_t dimension_;
    size_t size_;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    root = nullptr;
    dimension_ = N;
    size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    deleteTree(root);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::deleteTree(KDTreeNode<N, ElemType>* r)
{
    if (!r) return;
    deleteTree(r->nodes[0]);
    deleteTree(r->nodes[1]);
    size_--;
    delete r;
}

template <size_t N, typename ElemType>
KDTreeNode<N, ElemType>* KDTree<N, ElemType>::copy(KDTreeNode<N, ElemType>* r)
{
    if (!r) return r;
    KDTreeNode<N, ElemType>* temp = new KDTreeNode<N, ElemType>(r->point, r->value);
    temp->nodes[0] = copy(r->nodes[0]);
    temp->nodes[1] = copy(r->nodes[1]);
    return temp;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) { // arbol que recibe referencia de otro subarbol
    root = copy(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) { // operador = para inicializar en arbol

    if (this != &rhs)
    {
        deleteTree(root);
        root = copy(rhs.root);
        dimension_ = rhs.dimension_;
        size_ = rhs.size_;
    }
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {

    return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    // TODO(me): Fill this in.
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {

    if (size_) return false;
    return true;
}

// find, No pudo usarse porque tuve problemas con el const.
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(const Point<N>& pt, KDTreeNode<N, ElemType>**& p) { // devuelve el puntero a puntero en la posicion del nodo

    int cont = 0;
    for (p = &root; (*p) && (*p)->point != pt; p = &((*p)->nodes[pt[cont % N] > (*p)->point[cont % N]]), cont++); // Si la coordenada actual es menor al punto, va a la derecha
    return *p != nullptr; // si (*p) lo encontró,  devuelve 1
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {

    KDTreeNode<N, ElemType>** p;
    int cont = 0;
    for (p = &root; (*p) && (*p)->point != pt; p = &((*p)->nodes[pt[cont % N] > (*p)->point[cont % N]]), cont++); // Si la coordenada actual es menor al punto, va a la derecha
    return *p != nullptr; // si (*p) lo encontró,  devuelve 1
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {

    KDTreeNode<N, ElemType>** p;
    int cont = 0;
    for (p = &root; (*p) && (*p)->point != pt; p = &((*p)->nodes[pt[cont % N] > (*p)->point[cont % N]]), cont++);
    if (!(*p))
    {
        *p = new KDTreeNode<N, ElemType>(pt, value); // Si no está el punto, lo añade y size aumenta en uno.
        size_++;
    }
    else // Si el punto ya existe en el KD-Tree, su valor es sobreescrito
    {
        (*p)->value = value;
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    KDTreeNode<N, ElemType>** p;
    int cont = 0;
    for (p = &root; (*p) && (*p)->point != pt; p = &((*p)->nodes[pt[cont % N] > (*p)->point[cont % N]]), cont++);
    if (!(*p)) // Si el punto no existe, se crea un nodo con el valor por defecto de ElemType
    {
        ElemType element = ElemType(); // Valor por defecto del tipo ElemType
        *p = new KDTreeNode<N, ElemType>(pt, element); // Si no está el punto, lo añade y size aumenta en uno.
        size_++;
    }
    return (*p)->value; //Retorna una referencia al valor asociado
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    KDTreeNode<N, ElemType>** p;
    int cont = 0;
    for (p = &root; (*p) && (*p)->point != pt; p = &((*p)->nodes[pt[cont % N] > (*p)->point[cont % N]]), cont++);
    if (!(*p)) // Si no encuentra el punto, lanza una excepción de fuera de rango
    {
        throw out_of_range(" Estas fuera del rango. ");
    }
    else
        return (*p)->value; // Retorna una referencia al valor asociado
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    KDTreeNode<N, ElemType>** p;
    int cont = 0;
    for (p = &root; (*p) && (*p)->point != pt; p = &((*p)->nodes[pt[cont % N] > (*p)->point[cont % N]]), cont++);
    if (!(*p)) // Si no encuentra el punto, lanza una excepción de fuera de rango
    {
        throw out_of_range(" Estas fuera del rango. ");
    }
    else
        return (*p)->value; // Retorna una referencia al valor asociado
}


template <size_t N, typename ElemType>
void KDTree<N, ElemType>::neighbors(const Point<N> key, KDTreeNode<N, ElemType>* current_node, vector<pair<double, KDTreeNode<N, ElemType>*>>& nearest_neighbors_candidates, int depth, int k) const
{
    if (!current_node)
    {
        return;
    }
    nearest_neighbors_candidates.push_back(make_pair(distance(current_node->point, key), current_node)); // se almacenan los vecinos cercanos con su distancia
    std::sort(nearest_neighbors_candidates.begin(), nearest_neighbors_candidates.end());
    if (nearest_neighbors_candidates.size() > k) nearest_neighbors_candidates.pop_back(); // guarda y actualiza los k vecinos constantemente

    int axis = depth % dimension_;
    bool right = false;
    if (key[axis] < current_node->point[axis])
    {
        right = true;
        neighbors(key, current_node->nodes[0], nearest_neighbors_candidates, ++depth, k);
    }
    else
    {
        right = false;
        neighbors(key, current_node->nodes[1], nearest_neighbors_candidates, ++depth, k);
    }

    if (fabs(current_node->point[axis] - key[axis]) < nearest_neighbors_candidates[0].first)
    {
        if (right)
        {
            neighbors(key, current_node->nodes[0], nearest_neighbors_candidates, ++depth, k);
        }
        else
        {
            neighbors(key, current_node->nodes[1], nearest_neighbors_candidates, ++depth, k);
        }
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const {

    vector<ElemType> values = knn_query(key, k); // recibe los k vecinos mas cercanos
    vector<pair<int, ElemType>> pair_values;
    typename vector<ElemType>::iterator it;
    
    for (it = values.begin(); it != values.end(); it++) {
        pair_values.push_back(make_pair(count(values.begin(), values.end(), *it), *it));
    }
    sort(pair_values.begin(), pair_values.end());
    reverse(pair_values.begin(), pair_values.end());

    ElemType new_element = pair_values[0].second;


    return new_element;
}

template <size_t N, typename ElemType>
vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key, size_t k) const {

    vector<pair<double, KDTreeNode<N, ElemType>*>> pares;
    vector<ElemType> values;
    neighbors(key, root, pares, 0, k);// pares regresa lleno de las distancias

    typename vector<pair<double, KDTreeNode<N, ElemType>*>>::iterator it;
    for (it = pares.begin(); it != pares.end(); ++it)
        values.push_back(((*it).second)->value);

    return values;
}


#endif  // SRC_KDTREE_HPP_
