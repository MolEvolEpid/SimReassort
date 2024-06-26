// -*- C++ -*-
// GENEALOGY class

#ifndef _GENEALOGY_H_
#define _GENEALOGY_H_

#include <string>
#include <cstring>
#include <utility>

#include "internal.h"
#include "nodeseq.h"
#include "inventory.h"

static const size_t MEMORY_MAX = (1<<28); // 256MB

//! Encodes a genealogy.

//! A genealogy consists of a sequence of nodes
//! and the current time.
template <size_t ndeme = 1>
class genealogy_t : public nodeseq_t {

private:
  
  // GENEALOGY member data:
  // - a counter of serial numbers
  // - an initial time
  // - the current time
  // - a sequence of nodes
  
  //! The next unique name.
  name_t _unique;      

// private:
  //! The initial time.
  slate_t _t0;
  //! The current time.
  slate_t _time;
  //! Has the genealogy been obscured?
  //! This is set to false by default; 'obscure()' sets to true.
  bool _obscured;
  //! Has the genealogy been hided from the reassortment evenst?
  //! This is set to false by default; 'hide()' sets to true.
  bool _hided;
  
private:

  //! get the next unique name
  name_t unique (void) {
    name_t u = _unique;
    _unique++;
    return u;
  };
  //! clean up
  void clean (void) {
    _unique = 0;
    _t0 = _time = R_NaReal;
    _obscured = false;
    _hided = false;
  };

public:
  // SERIALIZATION
  //! size of serialized binary form
  size_t bytesize (void) const {
    return 3*sizeof(name_t) +
      2*sizeof(slate_t) + nodeseq_t::bytesize();
  };
  //! binary serialization
  friend raw_t* operator>> (const genealogy_t& G, raw_t* o) {
    name_t A[3]; A[0] = G._unique; A[1] = name_t(G._obscured); A[2] = name_t(G._hided);
    slate_t B[2]; B[0] = G._t0; B[1] = G._time;
    memcpy(o,A,sizeof(A)); o += sizeof(A);
    memcpy(o,B,sizeof(B)); o += sizeof(B);
    return reinterpret_cast<const nodeseq_t&>(G) >> o;
  };
  //! binary deserialization
  friend raw_t* operator>> (raw_t* o, genealogy_t& G) {
    G.clean();
    name_t A[3];
    slate_t B[2];
    memcpy(A,o,sizeof(A)); o += sizeof(A);
    memcpy(B,o,sizeof(B)); o += sizeof(B);
    G._unique = A[0]; G._obscured = bool(A[1]); G._hided = bool(A[2]);
    G._t0 = B[0]; G._time = B[1];
    return o >> reinterpret_cast<nodeseq_t&>(G);
  };

public:
  // CONSTRUCTORS
  //! basic constructor for genealogy class
  //!  t0 = initial time
  genealogy_t (double t0 = 0) {
    clean();
    _time = _t0 = slate_t(t0);
    _obscured = false;
    _hided = false;
  };
  //! constructor from serialized binary form
  genealogy_t (raw_t *o) {
    o >> *this;
  };
  //! copy constructor
  genealogy_t (const genealogy_t& G) {
    raw_t *o = new raw_t[G.bytesize()];
    G >> o >> *this;
    delete[] o;
  };
  //! copy assignment operator
  genealogy_t & operator= (const genealogy_t& G) {
    clean();
    raw_t *o = new raw_t[G.bytesize()];
    G >> o >> *this;
    delete[] o;
    return *this;
  };
  //! move constructor
  genealogy_t (genealogy_t&& G) {
    nodeseq_t::operator=(reinterpret_cast<const nodeseq_t&>(G));
    _unique = G._unique;
    _t0 = G._t0;
    _time = G._time;
    _obscured = G._obscured;
    _hided = G._hided;
  };
  //! move assignment operator
  genealogy_t& operator= (genealogy_t&& G) {
    clean();
    nodeseq_t::operator=(reinterpret_cast<const nodeseq_t&>(G));
    _unique = G._unique;
    _t0 = G._t0;
    _time = G._time;
    _obscured = G._obscured;
    _hided = G._hided;
    return *this;
  };
  //! destructor
  ~genealogy_t (void) {
    clean();
  };
  
  //! view/set current time.
  slate_t& time (void) {
    return _time;
  };
  //! view current time.
  slate_t time (void) const {
    return _time;
  };
  //! get zero time.
  slate_t timezero (void) const {
    return _t0;
  };
  //! update unique name
  void update_uniq (void) {
    _unique++;
  }

public:

  void lineage_count (double *t, int *ell) const {
    slate_t tcur = *t = timezero();
    size_t nd = (obscured()) ? 1 : ndeme;
    for (size_t j = 0; j < nd; j++) ell[j] = 0;
    for (node_it i = begin(); i != end(); i++) {
      if (tcur < (*i)->slate) {
        t++; ell += nd;
        *t = tcur = (*i)->slate;
        for (size_t j = 0; j < nd; j++) ell[j] = (ell-nd)[j];
      }
      (*i)->lineage_incr(ell);
    }
    t++; ell += nd;
    *t = time();
    for (size_t j = 0; j < nd; j++) ell[j] = 0;
  };

  SEXP lineage_count (void) const {
    SEXP t, ell, out, outn;
    int nt = ntime(timezero())+1;
    int nl = (obscured()) ? nt : ndeme*nt;
    PROTECT(t = NEW_NUMERIC(nt));
    PROTECT(ell = NEW_INTEGER(nl));
    PROTECT(out = NEW_LIST(2));
    PROTECT(outn = NEW_CHARACTER(2));
    set_list_elem(out,outn,t,"time",0);
    set_list_elem(out,outn,ell,"count",1);
    SET_NAMES(out,outn);
    lineage_count(REAL(t),INTEGER(ell));
    UNPROTECT(4);
    return out;
  };

public:
  
  //! R list description
  SEXP structure (void) const {
    SEXP O, On, T0, Time, Nodes;
    PROTECT(O = NEW_LIST(3));
    PROTECT(On = NEW_CHARACTER(3));
    PROTECT(Time = NEW_NUMERIC(1));
    *REAL(Time) = double(time());
    PROTECT(T0 = NEW_NUMERIC(1));
    *REAL(T0) = double(timezero());
    PROTECT(Nodes = nodeseq_t::structure());
    set_list_elem(O,On,Time,"time",0);
    set_list_elem(O,On,T0,"t0",1);
    set_list_elem(O,On,Nodes,"nodes",2);
    SET_NAMES(O,On);
    UNPROTECT(5);
    return O;
  };
  
public:

  //! human-readable info
  std::string describe (void) const {
    std::string o = "time = " + std::to_string(time()) + "\n"
      + "t0 = " + std::to_string(timezero()) + "\n"
      + nodeseq_t::describe();
    return o;
  };

public:
  
  //! machine-readable info
  virtual std::string yaml (std::string tab = "") const {
    std::string o;
    std::string t = tab + "  ";
    o = tab + "t0: " + std::to_string(timezero()) + "\n"
      + tab + "time: " + std::to_string(time()) + "\n"
      + tab + "nodes:\n" + nodeseq_t::yaml(tab);
    return o;
  };

public:

  //! put genealogy at current time into Newick format.
  std::string newick (bool compact = true) const {
    return nodeseq_t::newick(time(),compact);
  };

public:

  //! check the validity of the genealogy.
  void valid (void) const {};

  bool check_genealogy_size (size_t grace = 0) const {
    static size_t maxq = MEMORY_MAX/(sizeof(node_t)+2*sizeof(ball_t));
    bool ok = true;
    if (size() > maxq+grace) {
      err("maximum genealogy size exceeded!");
    } else if (size() > maxq) {
      ok = false;
    }
    return ok;
  };

private:

  node_t* make_node (color_t col, name_t d = 0) {
    check_genealogy_size(0);
    name_t u = unique();
    node_t *p = new node_t(u,_time,d);
    ball_t *g = new ball_t(p,u,green,d);
    ball_t *b = new ball_t(p,u,col,d);
    p->green_ball() = g;
    p->insert(g);
    p->insert(b);
    return p;
  };

  void destroy_node (node_t *p) {
    if (!p->holds_own())
      err("cannot destroy a node that does not hold its own green ball."); // # nocov
    if (p->size() != 2)
      err("cannot destroy a node with more than 2 balls."); // # nocov
    remove(p);
    delete p;
  };
  
  //! swap balls a and b, wherever they lie
  void swap (ball_t *a, ball_t *b) {
    node_t *p = a->holder();
    node_t *q = b->holder();
    if (p != q) {
      p->erase(a); q->insert(a); a->holder() = q;
      q->erase(b); p->insert(b); b->holder() = p;
    }
  };

  //! add node p; take as parent the node holding ball a.
  //! the deme of p is changed to match that of a
  void add (node_t *p, ball_t *a) {
    swap(a,p->green_ball());
    p->deme = a->deme();
    push_back(p);
  };
  
  // //! insert a node q after node p;
  // void insert_after (node_t *p, node_t *q) {
  //   for (node_it i = begin(); i != end(); i++) {
  //     if ((*i)->slate > p->slate) {
  //       insert(i,q);
  //       break;
  //     }
  //   }
  // };
  
  //! drop the node holding black ball a.
  void drop (ball_t *a) {
    if (!a->is(black))
      err("in '%s': inconceivable! (color: %s)",__func__,colores[a->color]); // #nocov
    node_t *p = a->holder();
    if (p->size() > 2) {   // pocket is large: we simply drop the ball
      p->erase(a);
      delete a;
    } else {  // pocket contains two: action depends on the other ball
      ball_t *b = p->other(a);
      switch (b->color) {
      case blue:                // change black ball for red ball
        a->color = red;
        break;
      case purple: case darkorange: case lightorange: // swap black ball for green ball, delete node
        a->deme() = p->deme;
        swap(a,p->green_ball());
        destroy_node(p);
        drop(a);                // recursively pursue dropping ball a
        break;
      case black: 
        if (b->deme() == p->deme) { // swap other for green, delete node
          swap(b,p->green_ball());
          destroy_node(p);
        } else {              // deme changes here, insert purple ball
          a->color = purple;
        }
        break;
      case green:
        if (b->owner()->deme == p->deme) { // swap other for green, delete node
          swap(b,p->green_ball());
          destroy_node(p);
        } else {
          a->color=purple;
        }
        break;
      case red: case grey:                             // #nocov
        err("in '%s': inconceivable error.",__func__); // #nocov
        break;
      }
    }
  };

public:
  //! birth into deme d 
  ball_t* birth (ball_t* a, slate_t t, name_t d = 0) {
    if ((!a->is(black)))
      err("in '%s': inconceivable! (1st color: %s)",__func__,colores[a->color]); // #nocov
    time() = t;
    node_t *p = make_node(black,d);
    ball_t *b = p->last_ball();
    p->slate = time();
    add(p,a);
    return b;           
  };
  //! birth of second or subsequent sibling into deme d
  ball_t* birth (node_t* p, name_t d = 0) {
    ball_t *b = new ball_t(p,unique(),black,d);
    p->insert(b);
    return b;
  };
  //! death
  void death (ball_t *a, slate_t t) {
    if ((!a->is(black)))
      err("in '%s': inconceivable! (1st color: %s)",__func__,colores[a->color]); // #nocov
    time() = t;
    drop(a);
  };
  //! graft a new lineage into deme d
  ball_t* graft (slate_t t, name_t d = 0) {
    time() = t;
    node_t *p = make_node(black,d);
    ball_t *b = p->last_ball();
    p->slate = timezero();
    push_front(p);
    return b;
  };
  //! insert a sample node
  void sample (ball_t* a, slate_t t, bool cont = false) {
    if ((!a->is(black)))
      err("in '%s': inconceivable! (1st color: %s)",__func__,colores[a->color]); // #nocov
    time() = t;
    node_t *p = make_node(blue,a->deme());
    p->slate = time();
    add(p,a);
    if (cont) a->color = red;
  };
  //! movement into deme d
  ball_t* migrate (ball_t* a, slate_t t, name_t d = 0) {
    time() = t;
    node_t *p = make_node(purple,a->deme());
    p->slate = time();
    add(p,a);
    a->deme() = d;
    return a;
  };
  //! reassort
  void reassort (ball_t *a, ball_t *b, slate_t t) {
    if ((!a->is(black)))
      err("in '%s': inconceivable! (1st color: %s)",__func__,colores[a->color]); // #nocov
    if ((!b->is(black)))
      err("in '%s': inconceivable! (1st color: %s)",__func__,colores[b->color]); // #nocov
    time() = t;
    node_t *p = a->holder();
    node_t *q = make_node(black,p->deme);
    ball_t *d = q->other(q->green_ball());
    swap(a,d);
    add(q,b);
    q->slate = t;
    //! only consider binary tree
    if (p->size() > 2) {
      err("Non-binary tree!");
    } else {
      drop(d);
      node_t *x = make_node(darkorange, a->deme());   // x: x + orange
      add(x,a);
      x->slate = t;
    }
  };
  //! indicate a reassortment event at that time
  //! in those un-reassorted segment trees
  void reassort_notice (ball_t *a, slate_t t) {
    if ((!a->is(black)))
      err("in '%s': inconceivable! (1st color: %s)",__func__,colores[a->color]); // #nocov
    time() = t;
    node_t *x = make_node(lightorange,a->deme());
    x->slate = t;
    add(x,a);
  };
  //! set up for extraction of black balls
  //! (see 'inventory.h')
  std::pair<node_it, node_it> extant (void) const {
    return std::pair<node_it,node_it>(cbegin(),cend());
  };

  //! prune the tree (drop all black balls)
  genealogy_t& prune (void) {
    pocket_t *blacks = colored(black);
    while (!blacks->empty()) {
      ball_t *b = *(blacks->begin());
      blacks->erase(b);
      drop(b);
    }
    delete blacks;
    return *this;
  };

  //! drop all purple balls
  //! and erase all deme information
  genealogy_t& obscure (void) {
    pocket_t *purples = colored(purple);
    while (!purples->empty()) {
      ball_t *a = *(purples->begin());
      node_t *p = a->holder();
      ball_t *b = p->other(a);
      swap(b,p->green_ball());
      destroy_node(p);
      purples->erase(a);
    }
    delete purples;
    pocket_t *blacks = colored(black);
    while (!blacks->empty()) {
      ball_t *a = *(blacks->begin());
      a->deme() = 0;
      blacks->erase(a);
    }
    delete blacks;
    for (node_it i = begin(); i != end(); i++)
      (*i)->deme = 0;
    _obscured = true;
    return *this;
  };

  bool obscured (void) const {
    return _obscured;
  };
  
  //! drop all orange balls
  genealogy_t& hide (void) {
    pocket_t *lightoranges = colored(lightorange);
    while (!lightoranges->empty()) {
      ball_t *a = *(lightoranges->begin());
      node_t *p = a->holder();
      ball_t *b = p->other(a);
      swap(b,p->green_ball());
      destroy_node(p);
      lightoranges->erase(a);
    }
    delete lightoranges;
    pocket_t *darkoranges = colored(darkorange);
    while (!darkoranges->empty()) {
      ball_t *a = *(darkoranges->begin());
      node_t *p = a->holder();
      ball_t *b = p->other(a);
      swap(b,p->green_ball());
      destroy_node(p);
      darkoranges->erase(a);
    }
    delete darkoranges;
    _hided = true;
    return *this;
  };
  
  bool hide (void) const {
    return _hided;
  };
  
  //! truncate the genealogy by removing nodes
  //! with times later than tnew
  //! NB: this destroys the genealogy inasmuch
  //! as the state is no longer correct.
  void truncate (slate_t tnew) {
    if (!empty()) {
      node_t *n = back();
      while (!empty() && n->slate > tnew) {
        if (n->holds(black)) {
          ball_t *b = n->last_ball(); // must be black!
          if (!b->is(black)) err("in '%s': inconceivable!",__func__); // #nocov
          drop(b);
        } else if (n->holds(red)) {
          ball_t *b = n->last_ball(); // must be red!
          if (!b->is(red)) err("in '%s': inconceivable!",__func__); // #nocov
          b->color = black;
          swap(b,n->green_ball());
          destroy_node(n);
        } else {
          err("in '%s': inconceivable error.",__func__); // #nocov
        }
        n = back();
      }
      time() = tnew;
    }
  };
};

#endif
