#ifndef H_BIT_MATRIX
#define H_BIT_MATRIX

#ifdef BITMATRIX_UNIT_TEST
#define BITMATRIX_CHECK(x,y,w,h) if (x>=w || y>=h) throw "Index out of bounds"
#else
#define BITMATRIX_CHECK(x,y,w,h)
#endif 

class BitMatrixRow
{
  typedef unsigned long long base;
  enum { BITS = (sizeof(base)*8) };
  typedef std::vector<base> Row;
  Row m_Row;
  friend class BitMatrix;

  static unsigned base_size(unsigned bits)
  {
    return (bits+BITS-1)/BITS;
  }

  base end_mask(unsigned width)
  {
    width=width&(BITS-1);
    if (width==0) return ~base(0);
    return (base(1)<<width)-1;
  }

  const base& get_word(unsigned x, unsigned& offset) const
  {
    unsigned word_index=x/BITS;
    offset=x&(BITS-1);
    return m_Row[word_index];
  }

  base& get_word(unsigned x, unsigned& offset)
  {
    unsigned word_index=x/BITS;
    offset=x&(BITS-1);
    return m_Row[word_index];
  }

  static base shift(unsigned offset)
  {
    base res=1;
    return res<<=offset;
  }
public:
  BitMatrixRow(unsigned width=0)
  {
    if (width>0)
      resize(width);
  }

  bool any() const
  {
    for(Row::const_iterator it=m_Row.begin();it!=m_Row.end();++it)
      if (*it) return true;
    return false;
  }

  void set()
  {
    std::fill(m_Row.begin(),m_Row.end(),~base(0));
  }

  void reset()
  {
    std::fill(m_Row.begin(),m_Row.end(),0);
  }

  bool test(unsigned x) const
  {
    unsigned offset;
    const base& word=get_word(x,offset);
    return ((word & shift(offset)) != 0);
  }

  void set(unsigned x, bool value=true)
  {
    unsigned offset;
    base& word=get_word(x,offset);
    if (value) word |= shift(offset);
    else       word &= ~shift(offset);
  }

  void reset(unsigned x)
  {
    unsigned offset;
    base& word=get_word(x,offset);
    word &= ~shift(offset);
  }

  void toggle(unsigned x)
  {
    unsigned offset;
    base& word=get_word(x,offset);
    base s=shift(offset);
    if ((word & s) != 0)
      word&=~s;
    else
      word|=s;
  }

  void resize(unsigned width)
  {
    m_Row.resize(base_size(width),0);
    m_Row.back() &= end_mask(width);
  }

  BitMatrixRow& operator &= (const BitMatrixRow& row)
  {
    Row::const_iterator i2=row.m_Row.begin(),e2=row.m_Row.end();
    for(Row::iterator i1=m_Row.begin();i1!=m_Row.end() && i2!=e2;++i1,++i2)
      *i1 &= *i2;
    return *this;
  }

  BitMatrixRow& operator ^= (const BitMatrixRow& row)
  {
    Row::const_iterator i2=row.m_Row.begin(),e2=row.m_Row.end();
    for(Row::iterator i1=m_Row.begin();i1!=m_Row.end() && i2!=e2;++i1,++i2)
      *i1 ^= *i2;
    return *this;
  }

  bool operator== (const BitMatrixRow& row) const
  {
    if (m_Row.size() != row.m_Row.size()) return false;
    return std::equal(m_Row.begin(),m_Row.end(),row.m_Row.begin());
  }

  bool operator!= (const BitMatrixRow& row) const
  {
    return !(*this == row);
  }

  const base& operator[] (unsigned idx) const { return m_Row[idx]; }
        base& operator[] (unsigned idx)       { return m_Row[idx]; }
};

inline bool operator== (const BitMatrixRow& a, const BitMatrixRow& b)
{

}

class BitMatrix
{
  typedef BitMatrixRow Row;
  typedef std::vector<Row> Matrix;
  typedef Row::base base;
  Matrix   m_Matrix;
  unsigned m_Width,m_Height;

public:
  BitMatrix(unsigned width=0, unsigned height=0)
  : m_Width(0)
  , m_Height(0)
  {
    resize(width,height);
  }

  unsigned size() const { return m_Matrix.size(); }

  typedef Matrix::const_iterator const_iterator;
  const_iterator begin() const { return m_Matrix.begin(); }
  const_iterator end()   const { return m_Matrix.end(); }

  BitMatrixRow& operator[] (unsigned y)
  {
    return m_Matrix[y];
  }

  const BitMatrixRow& operator[] (unsigned y) const
  {
    return m_Matrix[y];
  }

  unsigned get_width()  const { return m_Width; }
  unsigned get_height() const { return m_Height; }

  void resize(unsigned width, unsigned height)
  {
    if (width>0 && height>0)
    {
      for(Matrix::iterator it=m_Matrix.begin();it!=m_Matrix.end();++it)
        it->resize(width);
      m_Width=width;
      m_Height=height;
      m_Matrix.resize(height,Row(width));
    }
  }

  bool test(unsigned x, unsigned y) const
  {
    BITMATRIX_CHECK(x,y,m_Width,m_Height);
    return m_Matrix[y].test(x);
  }

  void set(unsigned x, unsigned y, bool value=true)
  {
    BITMATRIX_CHECK(x,y,m_Width,m_Height);
    m_Matrix[y].set(x,value);
  }

  void reset(unsigned x, unsigned y)
  {
    BITMATRIX_CHECK(x,y,m_Width,m_Height);
    m_Matrix[y].reset(x);
  }

  void toggle(unsigned x, unsigned y)
  {
    BITMATRIX_CHECK(x,y,m_Width,m_Height);
    m_Matrix[y].toggle(x);
  }

};

#endif // H_BIT_MATRIX

