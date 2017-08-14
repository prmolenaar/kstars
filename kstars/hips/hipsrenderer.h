/*
  Copyright (C) 2015-2017, Pavel Mraz

  Copyright (C) 2017, Jasem Mutlaq

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef HIPSRENDERER_H
#define HIPSRENDERER_H

#include "cmapview.h"
#include "cskpainter.h"
#include "hipsmanager.h"
#include "healpix.h"

class HiPSRenderer : public QObject
{
  Q_OBJECT
public:
  explicit HiPSRenderer();
  void render(mapView_t *view, CSkPainter *painter, QImage *pDest);
  void renderRec(bool allsky, int level, int pix, CSkPainter *painter, QImage *pDest);
  bool renderPix(bool allsky, int level, int pix, CSkPainter *painter, QImage *pDest);
  void setParam(const hipsParams_t &param);
  hipsParams_t *getParam();

  HiPSManager *manager();

signals:

public slots:

private:
  HiPSManager m_manager;
  int         m_blocks;
  int         m_rendered;
  int         m_size;
  QSet <int>  m_renderedMap;
  HEALPix     m_HEALpix;
};

extern HiPSRenderer *g_hipsRenderer;

#endif // HIPSRENDERER_H
