"""Script 25: interactive three volume viewer

Purpose
-------
Build the Three.js 3D viewer with embedded volume data and CDN-loaded modules.

Method overview
---------------
1. Read the smoothed dust cube in batches, determine the Cartesian grid and encode the
   configured density range as an unsigned 8-bit volume texture.
2. Build the fitted superbubble line geometry and Grade A/B bubble annotations,
   preferring the compact catalogue for display when it is available.
3. Embed the texture and annotations in the HTML viewer. Three.js modules load from a
   CDN; a separate binary payload is written only when --output-bin is supplied.

Main inputs
-----------
- ../results/intermediate_output/3d_dust_map_products/smoothed_3d_dust_cube.parquet
- ../data/Bubbles.csv
- ../results/superbubble_final_fit_parameters.csv
- ../results/Open_superbubbles.csv (optional simplified display table; the viewer falls back to the full table if it is absent)

Main outputs
------------
- ../results/interactive_3d_viewer.html
- Optional binary volume export with --output-bin.

Figure/table role
-----------------
../results/interactive_3d_viewer.html; local source for the OSB 3D visualization.

Runtime and data notes
----------------------
Reference runtime: 7.05 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
This script generates ../results/interactive_3d_viewer.html; the companion online platform is https://nadc.china-vo.org/data/dustmaps/OSBs.

Reading the code
----------------
Start with the path and scientific settings below, then follow main() at
the end of the file. The method overview above describes the order of the
analysis steps; the input/output lists identify the upstream data and
products written by this stage. Relative data paths are resolved from code/.
"""

from __future__ import annotations
import time

import argparse
import base64
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
FINAL_DIR = Path("..") / "results"

DUST_PARQUET = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products" / "smoothed_3d_dust_cube.parquet"
BUBBLES_CSV = DATA_DIR / "Bubbles.csv"
FINAL_TABLE = FINAL_DIR / "superbubble_final_fit_parameters.csv"
OPEN_SUPERBUBBLES_CSV = FINAL_DIR / "Open_superbubbles.csv"

OUTPUT_HTML = FINAL_DIR / "interactive_3d_viewer.html"

DEFAULT_X_RANGE = (-3.0, 3.0)
DEFAULT_Y_RANGE = (-3.0, 3.0)
DEFAULT_Z_RANGE = (-1.0, 1.0)
DEFAULT_DUST_TRANSPARENT_MAX = 0.05
DEFAULT_DUST_MAX = 0.6
DEFAULT_RAY_STEPS = 96
DEFAULT_OPACITY_SCALE = 4.0
DEFAULT_BRIGHTNESS = 1.45
PARQUET_BATCH_SIZE = 1_000_000


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Three.js Dust Volume Viewer</title>
  <style>
    html, body {
      margin: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      background: #07080a;
      font-family: Arial, Helvetica, sans-serif;
    }
    #app { width: 100vw; height: 100vh; }
    .toolbar {
      position: fixed;
      top: 10px;
      left: 10px;
      z-index: 10;
      display: flex;
      flex-direction: column;
      gap: 7px;
      padding: 8px;
      border: 1px solid rgba(255,255,255,0.16);
      border-radius: 6px;
      color: #eef2f5;
      background: rgba(9, 11, 15, 0.84);
      backdrop-filter: blur(6px);
      font-size: 12px;
    }
    .button-row {
      display: flex;
      gap: 6px;
      align-items: center;
      flex-wrap: nowrap;
    }
    .toolbar button {
      height: 28px;
      padding: 0 10px;
      border: 1px solid rgba(255,255,255,0.20);
      border-radius: 5px;
      color: #eef2f5;
      background: #26313d;
      font-size: 12px;
      cursor: pointer;
    }
    .toolbar button.is-off {
      color: #aab3bd;
      background: #161b22;
      border-color: rgba(255,255,255,0.10);
    }
    .control {
      display: grid;
      grid-template-columns: 72px 132px 38px;
      gap: 8px;
      align-items: center;
    }
    .control label {
      color: rgba(238,242,245,0.78);
      white-space: nowrap;
    }
    .control input {
      width: 132px;
      accent-color: #8bb7e0;
    }
    .control output {
      min-width: 38px;
      text-align: right;
      color: rgba(238,242,245,0.88);
      font-variant-numeric: tabular-nums;
    }
    .hud {
      position: fixed;
      right: 12px;
      bottom: 12px;
      z-index: 10;
      color: rgba(238,242,245,0.78);
      font-size: 12px;
      line-height: 1.4;
      text-align: right;
      pointer-events: none;
    }
    .tooltip {
      position: fixed;
      z-index: 20;
      display: none;
      max-width: 220px;
      padding: 5px 7px;
      border: 1px solid rgba(255,255,255,0.18);
      border-radius: 5px;
      color: #eef2f5;
      background: rgba(9, 11, 15, 0.90);
      font-size: 12px;
      line-height: 1.25;
      pointer-events: none;
      white-space: nowrap;
    }
  </style>
</head>
<body>
  <div class="toolbar">
    <div class="button-row">
      <button type="button" id="dustToggle">Dust</button>
      <button type="button" id="resetView">Reset view</button>
      <button type="button" id="osbToggle" class="is-off">Open superbubbles</button>
      <button type="button" id="bubblesToggle" class="is-off">Bubbles</button>
      <button type="button" id="rwToggle" class="is-off">Radcliffe Wave</button>
    </div>
    <div class="control">
      <label for="opacityScale">Opacity</label>
      <input id="opacityScale" type="range" min="0.10" max="15.00" step="0.05" value="__OPACITY_SCALE__">
      <output id="opacityValue">__OPACITY_SCALE__</output>
    </div>
    <div class="control">
      <label for="raySteps">Steps</label>
      <input id="raySteps" type="range" min="32" max="512" step="16" value="__RAY_STEPS__">
      <output id="stepsValue">__RAY_STEPS__</output>
    </div>
    <div class="control">
      <label for="brightness">Brightness</label>
      <input id="brightness" type="range" min="0.50" max="3.00" step="0.05" value="__BRIGHTNESS__">
      <output id="brightnessValue">__BRIGHTNESS__</output>
    </div>
  </div>
  <div id="app"></div>
  <div class="hud" id="hud"></div>
  <div class="tooltip" id="tooltip"></div>

  <script type="importmap">
    {
      "imports": {
        "three": "https://unpkg.com/three@0.160.0/build/three.module.js",
        "three/addons/": "https://unpkg.com/three@0.160.0/examples/jsm/"
      }
    }
  </script>
  <script type="module">
    import * as THREE from 'three';
    import { TrackballControls } from 'three/addons/controls/TrackballControls.js';

    const metadata = __METADATA_JSON__;
    const volumeBase64 = '__VOLUME_BASE64__';

    function decodeBase64(data) {
      const binary = atob(data);
      const bytes = new Uint8Array(binary.length);
      for (let i = 0; i < binary.length; i += 1) bytes[i] = binary.charCodeAt(i);
      return bytes;
    }

    const container = document.getElementById('app');
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x07080a);

    const camera = new THREE.PerspectiveCamera(55, window.innerWidth / window.innerHeight, 0.01, 100);
    const defaultCameraPosition = new THREE.Vector3(0.05, -0.05, 8.2);
    const defaultTarget = new THREE.Vector3(0, 0, 0);
    camera.position.copy(defaultCameraPosition);

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
    renderer.setPixelRatio(1);
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.domElement.style.touchAction = 'none';
    container.appendChild(renderer.domElement);

    const controls = new TrackballControls(camera, renderer.domElement);
    controls.target.copy(defaultTarget);
    controls.rotateSpeed = 3.2;
    controls.zoomSpeed = 1.1;
    controls.panSpeed = 0.75;
    controls.staticMoving = false;
    controls.dynamicDampingFactor = 0.12;
    controls.update();

    const volumeData = decodeBase64(volumeBase64);
    const volumeTexture = new THREE.Data3DTexture(
      volumeData,
      metadata.grid.nx,
      metadata.grid.ny,
      metadata.grid.nz
    );
    volumeTexture.format = THREE.RedFormat;
    volumeTexture.type = THREE.UnsignedByteType;
    volumeTexture.minFilter = THREE.LinearFilter;
    volumeTexture.magFilter = THREE.LinearFilter;
    volumeTexture.wrapS = THREE.ClampToEdgeWrapping;
    volumeTexture.wrapT = THREE.ClampToEdgeWrapping;
    volumeTexture.wrapR = THREE.ClampToEdgeWrapping;
    volumeTexture.unpackAlignment = 1;
    volumeTexture.needsUpdate = true;

    const boxMin = new THREE.Vector3(metadata.bounds.xMin, metadata.bounds.yMin, metadata.bounds.zMin);
    const boxMax = new THREE.Vector3(metadata.bounds.xMax, metadata.bounds.yMax, metadata.bounds.zMax);
    const boxSize = new THREE.Vector3().subVectors(boxMax, boxMin);
    const boxCenter = new THREE.Vector3().addVectors(boxMin, boxMax).multiplyScalar(0.5);

    const volumeGeometry = new THREE.BoxGeometry(boxSize.x, boxSize.y, boxSize.z);
    volumeGeometry.translate(boxCenter.x, boxCenter.y, boxCenter.z);

    const material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3,
      side: THREE.BackSide,
      transparent: true,
      depthWrite: false,
      uniforms: {
        uData: { value: volumeTexture },
        uBoxMin: { value: boxMin },
        uBoxMax: { value: boxMax },
        uSteps: { value: metadata.render.raySteps },
        uOpacityScale: { value: metadata.render.opacityScale },
        uBrightness: { value: metadata.render.brightness }
      },
      vertexShader: `
        out vec3 vPosition;
        void main() {
          vPosition = position;
          gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
        }
      `,
      fragmentShader: `
        precision highp float;
        precision highp sampler3D;

        in vec3 vPosition;
        uniform sampler3D uData;
        uniform vec3 uBoxMin;
        uniform vec3 uBoxMax;
        uniform int uSteps;
        uniform float uOpacityScale;
        uniform float uBrightness;
        out vec4 outColor;

        vec2 intersectBox(vec3 origin, vec3 direction) {
          vec3 invDirection = 1.0 / direction;
          vec3 t0 = (uBoxMin - origin) * invDirection;
          vec3 t1 = (uBoxMax - origin) * invDirection;
          vec3 tMin = min(t0, t1);
          vec3 tMax = max(t0, t1);
          float nearT = max(max(tMin.x, tMin.y), tMin.z);
          float farT = min(min(tMax.x, tMax.y), tMax.z);
          return vec2(nearT, farT);
        }

        void main() {
          vec3 rayDirection = normalize(vPosition - cameraPosition);
          vec2 hit = intersectBox(cameraPosition, rayDirection);
          if (hit.x > hit.y) discard;

          float startT = max(hit.x, 0.0);
          float rayLength = hit.y - startT;
          if (rayLength <= 0.0) discard;

          int steps = clamp(uSteps, 8, 512);
          float delta = rayLength / float(steps);
          vec3 position = cameraPosition + rayDirection * startT;
          vec3 stepVector = rayDirection * delta;
          vec3 boxScale = uBoxMax - uBoxMin;
          vec4 accumulated = vec4(0.0);

          for (int i = 0; i < 512; i += 1) {
            if (i >= steps) break;

            vec3 texCoord = (position - uBoxMin) / boxScale;
            float density = texture(uData, texCoord).r;
            if (density > 0.0) {
              float sampleAlpha = 1.0 - exp(-density * uOpacityScale * delta);
              vec3 sampleColor = vec3(clamp(pow(density, 0.75) * uBrightness, 0.0, 1.0));
              accumulated.rgb += (1.0 - accumulated.a) * sampleAlpha * sampleColor;
              accumulated.a += (1.0 - accumulated.a) * sampleAlpha;
              if (accumulated.a > 0.97) break;
            }

            position += stepVector;
          }

          if (accumulated.a <= 0.002) discard;
          outColor = accumulated;
        }
      `
    });

    const volumeMesh = new THREE.Mesh(volumeGeometry, material);
    scene.add(volumeMesh);

    const edges = new THREE.LineSegments(
      new THREE.EdgesGeometry(volumeGeometry),
      new THREE.LineBasicMaterial({ color: 0x48627f, transparent: true, opacity: 0.34 })
    );
    scene.add(edges);

    const axes = new THREE.AxesHelper(1.0);
    axes.material.depthTest = false;
    axes.renderOrder = 2;
    scene.add(axes);

    const grid = new THREE.GridHelper(6, 12, 0x31506f, 0x223348);
    grid.rotation.x = Math.PI / 2;
    grid.material.opacity = 0.20;
    grid.material.transparent = true;
    scene.add(grid);

    const pickableObjects = [];

    const osbRoot = new THREE.Group();
    osbRoot.visible = false;
    scene.add(osbRoot);
    for (const group of metadata.openSuperbubbles || []) {
      for (const item of group.items || []) {
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(item.segments.flat(2), 3));
        const lineMaterial = new THREE.LineBasicMaterial({
          color: new THREE.Color(group.color),
          transparent: true,
          opacity: 0.95,
          depthWrite: false
        });
        const lines = new THREE.LineSegments(geometry, lineMaterial);
        lines.renderOrder = 3;
        lines.userData.tooltip = 'Open superbubble ' + item.id;
        osbRoot.add(lines);
        pickableObjects.push(lines);
      }
    }

    const rwRoot = new THREE.Group();
    rwRoot.visible = false;
    scene.add(rwRoot);
    const rwSlope = Math.tan(Math.PI / 3);
    const rwIntercept = 0.6;
    const rwXMin = Math.max(metadata.bounds.xMin, (metadata.bounds.yMin - rwIntercept) / rwSlope);
    const rwXMax = Math.min(metadata.bounds.xMax, (metadata.bounds.yMax - rwIntercept) / rwSlope);
    const rwYMin = rwSlope * rwXMin + rwIntercept;
    const rwYMax = rwSlope * rwXMax + rwIntercept;
    const rwGeometry = new THREE.BufferGeometry();
    rwGeometry.setAttribute('position', new THREE.Float32BufferAttribute([
      rwXMin, rwYMin, metadata.bounds.zMin,
      rwXMax, rwYMax, metadata.bounds.zMin,
      rwXMax, rwYMax, metadata.bounds.zMax,
      rwXMin, rwYMin, metadata.bounds.zMax
    ], 3));
    rwGeometry.setIndex([0, 1, 2, 0, 2, 3]);
    rwGeometry.computeVertexNormals();
    const rwMaterial = new THREE.MeshBasicMaterial({
      color: 0xffb000,
      transparent: true,
      opacity: 0.22,
      side: THREE.DoubleSide,
      depthWrite: false
    });
    const rwPlane = new THREE.Mesh(rwGeometry, rwMaterial);
    rwPlane.renderOrder = 2;
    rwRoot.add(rwPlane);
    const rwEdges = new THREE.LineSegments(
      new THREE.EdgesGeometry(rwGeometry),
      new THREE.LineBasicMaterial({ color: 0xffc857, transparent: true, opacity: 0.85, depthWrite: false })
    );
    rwRoot.add(rwEdges);

    const bubblesRoot = new THREE.Group();
    bubblesRoot.visible = false;
    scene.add(bubblesRoot);
    const bubbleGeometry = new THREE.SphereGeometry(1, 20, 10);
    const bubbleMaterial = new THREE.MeshBasicMaterial({
      color: 0xf6c343,
      transparent: true,
      opacity: 0.55,
      wireframe: true,
      depthWrite: false
    });
    const bubbles = metadata.bubbles || [];
    const bubbleMesh = new THREE.InstancedMesh(bubbleGeometry, bubbleMaterial, bubbles.length);
    bubbleMesh.frustumCulled = false;
    const bubbleObject = new THREE.Object3D();
    bubbles.forEach((bubble, index) => {
      const radius = Math.max(0.001, bubble.radius_kpc);
      bubbleObject.position.set(bubble.x, bubble.y, bubble.z);
      bubbleObject.scale.set(radius, radius, radius);
      bubbleObject.updateMatrix();
      bubbleMesh.setMatrixAt(index, bubbleObject.matrix);
    });
    bubbleMesh.instanceMatrix.needsUpdate = true;
    bubblesRoot.add(bubbleMesh);
    pickableObjects.push(bubbleMesh);

    const openSuperbubbleCount = (metadata.openSuperbubbles || [])
      .reduce((sum, group) => sum + (group.items || []).length, 0);

    document.getElementById('hud').innerHTML =
      'Three.js ray marching<br>' +
      metadata.grid.nx.toLocaleString() + ' x ' +
      metadata.grid.ny.toLocaleString() + ' x ' +
      metadata.grid.nz.toLocaleString() + ' volume texture<br>' +
      metadata.grid.totalVoxels.toLocaleString() + ' sampled cells<br>' +
      metadata.transfer.dustTransparentMax.toFixed(2) + '-' +
      metadata.transfer.dustMax.toFixed(2) + ' mag/kpc mapped to density<br>' +
      openSuperbubbleCount.toLocaleString() + ' Open superbubbles<br>' +
      bubbles.length.toLocaleString() + ' Grade A/B bubbles<br>' +
      'Radcliffe Wave plane: y = tan(60 deg) x + 0.6 kpc';

    const dustToggle = document.getElementById('dustToggle');
    dustToggle.addEventListener('click', () => {
      volumeMesh.visible = !volumeMesh.visible;
      dustToggle.classList.toggle('is-off', !volumeMesh.visible);
    });

    const osbToggle = document.getElementById('osbToggle');
    osbToggle.addEventListener('click', () => {
      osbRoot.visible = !osbRoot.visible;
      osbToggle.classList.toggle('is-off', !osbRoot.visible);
    });

    const bubblesToggle = document.getElementById('bubblesToggle');
    bubblesToggle.addEventListener('click', () => {
      bubblesRoot.visible = !bubblesRoot.visible;
      bubblesToggle.classList.toggle('is-off', !bubblesRoot.visible);
    });

    const rwToggle = document.getElementById('rwToggle');
    rwToggle.addEventListener('click', () => {
      rwRoot.visible = !rwRoot.visible;
      rwToggle.classList.toggle('is-off', !rwRoot.visible);
    });

    document.getElementById('resetView').addEventListener('click', () => {
      camera.position.copy(defaultCameraPosition);
      camera.up.set(0, 1, 0);
      controls.target.copy(defaultTarget);
      controls.update();
    });

    const opacityInput = document.getElementById('opacityScale');
    const opacityOutput = document.getElementById('opacityValue');
    opacityInput.addEventListener('input', () => {
      const value = parseFloat(opacityInput.value);
      material.uniforms.uOpacityScale.value = value;
      opacityOutput.value = value.toFixed(2);
    });

    const stepsInput = document.getElementById('raySteps');
    const stepsOutput = document.getElementById('stepsValue');
    stepsInput.addEventListener('input', () => {
      const value = parseInt(stepsInput.value, 10);
      material.uniforms.uSteps.value = value;
      stepsOutput.value = String(value);
    });

    const brightnessInput = document.getElementById('brightness');
    const brightnessOutput = document.getElementById('brightnessValue');
    brightnessInput.addEventListener('input', () => {
      const value = parseFloat(brightnessInput.value);
      material.uniforms.uBrightness.value = value;
      brightnessOutput.value = value.toFixed(2);
    });

    const tooltip = document.getElementById('tooltip');
    const pointer = new THREE.Vector2();
    const raycaster = new THREE.Raycaster();
    raycaster.params.Line.threshold = 0.035;

    function hideTooltip() {
      tooltip.style.display = 'none';
    }

    function showTooltip(event, text) {
      tooltip.textContent = text;
      tooltip.style.left = String(event.clientX + 12) + 'px';
      tooltip.style.top = String(event.clientY + 12) + 'px';
      tooltip.style.display = 'block';
    }

    renderer.domElement.addEventListener('pointermove', (event) => {
      const rect = renderer.domElement.getBoundingClientRect();
      pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
      pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
      raycaster.setFromCamera(pointer, camera);

      const intersections = raycaster.intersectObjects(pickableObjects, false);
      for (const hit of intersections) {
        if (hit.object === bubbleMesh && bubblesRoot.visible && Number.isInteger(hit.instanceId)) {
          const bubble = bubbles[hit.instanceId];
          if (bubble) {
            showTooltip(event, 'Bubble ' + bubble.id + '  Grade ' + bubble.grade);
            return;
          }
        }
        if (hit.object.userData.tooltip && hit.object.parent && hit.object.parent.visible) {
          showTooltip(event, hit.object.userData.tooltip);
          return;
        }
      }
      hideTooltip();
    });

    renderer.domElement.addEventListener('pointerleave', hideTooltip);

    window.addEventListener('resize', () => {
      camera.aspect = window.innerWidth / window.innerHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(window.innerWidth, window.innerHeight);
      if (typeof controls.handleResize === 'function') controls.handleResize();
      controls.update();
    });

    function animate() {
      requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene, camera);
    }
    animate();
  </script>
</body>
</html>
"""


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError('Missing input file.')


def finite_float(value, default=None):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if np.isfinite(number):
        return number
    return default


def rotate_xy(x, y, angle_deg):
    angle = np.deg2rad(angle_deg)
    ca = np.cos(angle)
    sa = np.sin(angle)
    return ca * x - sa * y, sa * x + ca * y


def add_polyline_segments(target, points):
    for left, right in zip(points[:-1], points[1:]):
        target.append(left)
        target.append(right)


def cap_segments(sb, n=72):
    bx = float(sb["base_x"])
    by = float(sb["base_y"])
    bz = float(sb["base_z"])
    a = float(sb["a"])
    b = float(sb["b"])
    c = float(sb["c_cap"])
    sign = float(sb["cap_sign"])
    angle_deg = float(sb["pa_deg"])
    theta = np.linspace(0, 2 * np.pi, n, endpoint=True)
    segments = []

    for frac in np.linspace(0.0, 1.0, 6):
        z_local = sign * c * frac
        scale = np.sqrt(max(0.0, 1.0 - frac * frac))
        x0 = a * scale * np.cos(theta)
        y0 = b * scale * np.sin(theta)
        xr, yr = rotate_xy(x0, y0, angle_deg)
        points = np.column_stack([bx + xr, by + yr, np.full_like(theta, bz + z_local)])
        add_polyline_segments(segments, points.tolist())

    z_frac = np.linspace(0.0, 1.0, n)
    for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        scale = np.sqrt(np.maximum(0.0, 1.0 - z_frac * z_frac))
        x0 = a * scale * np.cos(phi)
        y0 = b * scale * np.sin(phi)
        xr, yr = rotate_xy(x0, y0, angle_deg)
        points = np.column_stack([bx + xr, by + yr, bz + sign * c * z_frac])
        add_polyline_segments(segments, points.tolist())
    return segments


def cylinder_segments(sb, n=72):
    cx = float(sb["x"])
    cy = float(sb["y"])
    cz = float(sb["z"])
    a = float(sb["a"])
    b = float(sb["b"])
    h = float(sb["half_height"])
    angle_deg = float(sb["pa_deg"])
    theta = np.linspace(0, 2 * np.pi, n, endpoint=True)
    x0 = a * np.cos(theta)
    y0 = b * np.sin(theta)
    xr, yr = rotate_xy(x0, y0, angle_deg)
    segments = []
    add_polyline_segments(segments, np.column_stack([cx + xr, cy + yr, np.full_like(theta, cz - h)]).tolist())
    add_polyline_segments(segments, np.column_stack([cx + xr, cy + yr, np.full_like(theta, cz + h)]).tolist())
    for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        x_side, y_side = rotate_xy(a * np.cos(phi), b * np.sin(phi), angle_deg)
        segments.append([cx + x_side, cy + y_side, cz - h])
        segments.append([cx + x_side, cy + y_side, cz + h])
    return segments


def galactic_to_xyz(l_deg: np.ndarray, b_deg: np.ndarray, distance_kpc: np.ndarray):
    l_rad = np.deg2rad(l_deg)
    b_rad = np.deg2rad(b_deg)
    cos_b = np.cos(b_rad)
    x = distance_kpc * cos_b * np.cos(l_rad)
    y = distance_kpc * cos_b * np.sin(l_rad)
    z = distance_kpc * np.sin(b_rad)
    return x, y, z


def build_superbubble_annotations(final_table: Path, open_csv: Path):
    geometry_table = pd.read_csv(final_table)
    geometry_lookup = {str(row.get("id")): row for _, row in geometry_table.iterrows()}
    display_table = pd.read_csv(open_csv) if open_csv.exists() else geometry_table

    rows = []
    for _, display_row in display_table.iterrows():
        sb_id = str(display_row.get("id"))
        if sb_id not in geometry_lookup:
            continue
        row = geometry_lookup[sb_id]
        shape = str(row.get("shape", "")).lower()
        mark = int(finite_float(row.get("mark"), 0) or 0)
        sb_class = str(display_row.get("class", ""))
        sb_type = str(display_row.get("type", ""))

        if shape == "ellipsoid":
            base_x = finite_float(row.get("xy_plane_center_x_kpc"))
            base_y = finite_float(row.get("xy_plane_center_y_kpc"))
            base_z = finite_float(row.get("xy_plane_z_kpc"))
            a = finite_float(row.get("xy_plane_a_kpc"))
            b = finite_float(row.get("xy_plane_b_kpc"))
            angle = finite_float(row.get("xy_plane_angle_deg"), 0.0)
            center_z = finite_float(row.get("center_z_kpc"))
            c_full = finite_float(row.get("c_radius_kpc"))
            if None in (base_x, base_y, base_z, a, b, center_z, c_full):
                continue
            if mark == 1:
                z_apex = center_z - c_full
                cap_sign = -1
            elif mark == 2:
                z_apex = center_z + c_full
                cap_sign = 1
            else:
                continue
            rows.append(
                {
                    "id": sb_id,
                    "class": sb_class,
                    "type": sb_type,
                    "shape": "ellipsoid_cap",
                    "mark": mark,
                    "base_x": base_x,
                    "base_y": base_y,
                    "base_z": base_z,
                    "a": a,
                    "b": b,
                    "c_cap": abs(base_z - z_apex),
                    "cap_sign": cap_sign,
                    "z_apex": z_apex,
                    "pa_deg": angle,
                }
            )
            continue

        if shape == "cylinder":
            x = finite_float(row.get("center_x_kpc"))
            y = finite_float(row.get("center_y_kpc"))
            z = finite_float(row.get("center_z_kpc"), finite_float(row.get("xy_plane_z_kpc"), 0.0))
            a = finite_float(row.get("a_radius_kpc"), finite_float(row.get("xy_plane_a_kpc")))
            b = finite_float(row.get("b_radius_kpc"), finite_float(row.get("xy_plane_b_kpc")))
            angle = finite_float(row.get("angle_deg"), finite_float(row.get("xy_plane_angle_deg"), 0.0))
            if None in (x, y, z, a, b):
                continue
            rows.append(
                {
                    "id": sb_id,
                    "class": sb_class,
                    "type": sb_type,
                    "shape": "elliptical_cylinder",
                    "mark": mark,
                    "x": x,
                    "y": y,
                    "z": z,
                    "a": a,
                    "b": b,
                    "half_height": 0.025,
                    "pa_deg": angle,
                }
            )
    return rows


def build_bubble_annotations(bubbles_csv: Path):
    bubbles = pd.read_csv(bubbles_csv)
    grade = bubbles["Grade"].astype(str).str.strip().str.upper()
    bubbles = bubbles.loc[grade.isin(["A", "B"])].copy()
    distance = pd.to_numeric(bubbles["best_distance"], errors="coerce")
    l_deg = pd.to_numeric(bubbles["l"], errors="coerce")
    b_deg = pd.to_numeric(bubbles["b"], errors="coerce")
    ok = np.isfinite(l_deg) & np.isfinite(b_deg) & np.isfinite(distance)
    bubbles = bubbles.loc[ok].copy()
    x, y, z = galactic_to_xyz(l_deg.loc[ok].to_numpy(), b_deg.loc[ok].to_numpy(), distance.loc[ok].to_numpy())
    radius_kpc = pd.to_numeric(bubbles["physical_size"], errors="coerce").fillna(20.0).to_numpy() / 2000.0

    rows = []
    for i, row in bubbles.reset_index(drop=True).iterrows():
        rows.append(
            {
                "id": str(row.get("ID")),
                "grade": str(row.get("Grade", "")),
                "x": float(x[i]),
                "y": float(y[i]),
                "z": float(z[i]),
                "radius_kpc": float(radius_kpc[i]),
            }
        )
    return rows


def build_open_superbubble_lines(superbubbles):
    groups = {
        "northOpenCaps": {"label": "North-open caps", "color": "#1e88e5", "items": []},
        "southOpenCaps": {"label": "South-open caps", "color": "#e53935", "items": []},
        "ellipticalCylinders": {"label": "Elliptical cylinders", "color": "#2e7d32", "items": []},
    }
    for sb in superbubbles:
        if sb.get("shape") == "ellipsoid_cap":
            key = "northOpenCaps" if int(sb.get("mark", 0)) == 1 else "southOpenCaps"
            segments = cap_segments(sb)
        elif sb.get("shape") == "elliptical_cylinder":
            key = "ellipticalCylinders"
            segments = cylinder_segments(sb)
        else:
            continue
        groups[key]["items"].append(
            {
                "id": str(sb.get("id", "")),
                "shape": str(sb.get("shape", "")),
                "segments": segments,
            }
        )
    return [group for group in groups.values() if group["items"]]


def parquet_column_names(parquet_path: Path):
    dataset = ds.dataset(parquet_path, format="parquet")
    names = set(dataset.schema.names)
    if {"x", "y", "z", "dust"}.issubset(names):
        return "x", "y", "z", "dust"
    if {"X", "Y", "Z", "dust"}.issubset(names):
        return "X", "Y", "Z", "dust"
    if {"x", "y", "z", "dust_raw"}.issubset(names):
        return "x", "y", "z", "dust_raw"
    if {"X", "Y", "Z", "dust_raw"}.issubset(names):
        return "X", "Y", "Z", "dust_raw"
    raise ValueError(f"Could not identify dust parquet fields: {dataset.schema.names}")


def make_filter(x_col, y_col, z_col, x_range, y_range, z_range):
    return (
        (ds.field(x_col) >= float(x_range[0]))
        & (ds.field(x_col) <= float(x_range[1]))
        & (ds.field(y_col) >= float(y_range[0]))
        & (ds.field(y_col) <= float(y_range[1]))
        & (ds.field(z_col) >= float(z_range[0]))
        & (ds.field(z_col) <= float(z_range[1]))
    )


def collect_axes(dataset, columns, filter_expr, batch_size):
    axis_sets = {name: set() for name in columns}
    scanner = dataset.scanner(columns=list(columns), filter=filter_expr, batch_size=batch_size)
    for batch in scanner.to_batches():
        for name in columns:
            values = batch.column(name).to_numpy(zero_copy_only=False)
            axis_sets[name].update(np.unique(values).astype(float).tolist())
    axes = {name: np.array(sorted(values), dtype=np.float32) for name, values in axis_sets.items()}
    for name, axis in axes.items():
        if axis.size == 0:
            raise ValueError('Invalid input or missing required data.')
    return axes


def axis_step(axis):
    if axis.size < 2:
        return 1.0
    diffs = np.diff(axis.astype(float))
    return float(np.median(diffs))


def build_volume_from_parquet(
    parquet_path: Path,
    x_range,
    y_range,
    z_range,
    dust_transparent_max: float,
    dust_max: float,
    batch_size: int,
):
    ensure_file(parquet_path, "smoothed 3D dust cube parquet")
    if dust_transparent_max < 0:
        raise ValueError('Invalid input or missing required data.')
    if dust_max <= dust_transparent_max:
        raise ValueError('Invalid input or missing required data.')

    dataset = ds.dataset(parquet_path, format="parquet")
    x_col, y_col, z_col, dust_col = parquet_column_names(parquet_path)
    filter_expr = make_filter(x_col, y_col, z_col, x_range, y_range, z_range)

    print("Scanning dust grid axes...")
    axes = collect_axes(dataset, (x_col, y_col, z_col), filter_expr, batch_size)
    x_axis = axes[x_col]
    y_axis = axes[y_col]
    z_axis = axes[z_col]
    nx, ny, nz = len(x_axis), len(y_axis), len(z_axis)
    print(f"  grid dimensions: {nx} x {ny} x {nz} = {nx * ny * nz:,} cells")

    dx = axis_step(x_axis)
    dy = axis_step(y_axis)
    dz = axis_step(z_axis)
    volume = np.zeros((nx, ny, nz), dtype=np.uint8)

    print("Writing the uint8 volume texture...")
    scanner = dataset.scanner(columns=[x_col, y_col, z_col, dust_col], filter=filter_expr, batch_size=batch_size)
    for batch in scanner.to_batches():
        bx = batch.column(x_col).to_numpy(zero_copy_only=False).astype(np.float32)
        by = batch.column(y_col).to_numpy(zero_copy_only=False).astype(np.float32)
        bz = batch.column(z_col).to_numpy(zero_copy_only=False).astype(np.float32)
        dust = batch.column(dust_col).to_numpy(zero_copy_only=False).astype(np.float32)

        ix = np.rint((bx - x_axis[0]) / dx).astype(np.int64)
        iy = np.rint((by - y_axis[0]) / dy).astype(np.int64)
        iz = np.rint((bz - z_axis[0]) / dz).astype(np.int64)
        valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (iz >= 0) & (iz < nz)
        if not np.any(valid):
            continue
        clean = np.nan_to_num(dust[valid], nan=0.0, posinf=dust_max, neginf=0.0)
        normalized = np.clip((clean - dust_transparent_max) / (dust_max - dust_transparent_max), 0.0, 1.0)
        volume[ix[valid], iy[valid], iz[valid]] = np.rint(normalized * 255.0).astype(np.uint8)

    texture_order = np.transpose(volume, (2, 1, 0)).copy(order="C")
    metadata = {
        "source": str(parquet_path.as_posix()),
        "grid": {
            "nx": int(nx),
            "ny": int(ny),
            "nz": int(nz),
            "totalVoxels": int(volume.size),
        },
        "bounds": {
            "xMin": float(np.min(x_axis)),
            "xMax": float(np.max(x_axis)),
            "yMin": float(np.min(y_axis)),
            "yMax": float(np.max(y_axis)),
            "zMin": float(np.min(z_axis)),
            "zMax": float(np.max(z_axis)),
        },
        "transfer": {
            "dustTransparentMax": float(dust_transparent_max),
            "dustMax": float(dust_max),
            "encoding": "uint8 density = clip((dust - transparent_max) / (dust_max - transparent_max), 0, 1)",
        },
    }
    return texture_order, metadata


def write_viewer(texture_order: np.ndarray, metadata: dict, output_html: Path, output_bin: Path | None = None):
    output_html.parent.mkdir(parents=True, exist_ok=True)
    if output_bin is not None:
        output_bin.parent.mkdir(parents=True, exist_ok=True)
        output_bin.write_bytes(texture_order.tobytes())
    html = (
        HTML_TEMPLATE.replace("__METADATA_JSON__", json.dumps(metadata, separators=(",", ":")))
        .replace("__VOLUME_BASE64__", base64.b64encode(texture_order.tobytes()).decode("ascii"))
        .replace("__RAY_STEPS__", str(int(metadata["render"]["raySteps"])))
        .replace("__OPACITY_SCALE__", f"{float(metadata['render']['opacityScale']):.2f}")
        .replace("__BRIGHTNESS__", f"{float(metadata['render']['brightness']):.2f}")
    )
    output_html.write_text(html, encoding="utf-8")


def parse_args():
    """Build the command-line interface for Script 25."""
    parser = argparse.ArgumentParser(description='Run Script 25: interactive three volume viewer.')
    parser.add_argument("--dust-parquet", default=str(DUST_PARQUET))
    parser.add_argument("--final-table", default=str(FINAL_TABLE))
    parser.add_argument("--open-superbubbles", default=str(OPEN_SUPERBUBBLES_CSV))
    parser.add_argument("--bubbles", default=str(BUBBLES_CSV))
    parser.add_argument("--output-html", default=str(OUTPUT_HTML))
    parser.add_argument("--output-bin", default=None,
                        help="Optional separate binary volume export; the HTML already embeds the volume.")
    parser.add_argument("--x-min", type=float, default=DEFAULT_X_RANGE[0])
    parser.add_argument("--x-max", type=float, default=DEFAULT_X_RANGE[1])
    parser.add_argument("--y-min", type=float, default=DEFAULT_Y_RANGE[0])
    parser.add_argument("--y-max", type=float, default=DEFAULT_Y_RANGE[1])
    parser.add_argument("--z-min", type=float, default=DEFAULT_Z_RANGE[0])
    parser.add_argument("--z-max", type=float, default=DEFAULT_Z_RANGE[1])
    parser.add_argument("--dust-transparent-max", type=float, default=DEFAULT_DUST_TRANSPARENT_MAX)
    parser.add_argument("--dust-max", type=float, default=DEFAULT_DUST_MAX)
    parser.add_argument("--ray-steps", type=int, default=DEFAULT_RAY_STEPS)
    parser.add_argument("--opacity-scale", type=float, default=DEFAULT_OPACITY_SCALE)
    parser.add_argument("--brightness", type=float, default=DEFAULT_BRIGHTNESS)
    parser.add_argument("--batch-size", type=int, default=PARQUET_BATCH_SIZE)
    return parser.parse_args()


def main():
    """Run Script 25 from validated inputs to the documented outputs."""
    args = parse_args()
    ensure_file(Path(args.final_table), "final superbubble parameter table")
    ensure_file(Path(args.bubbles), "Wang et al. 2026 bubble catalog")

    texture_order, metadata = build_volume_from_parquet(
        Path(args.dust_parquet),
        (args.x_min, args.x_max),
        (args.y_min, args.y_max),
        (args.z_min, args.z_max),
        args.dust_transparent_max,
        args.dust_max,
        args.batch_size,
    )
    superbubbles = build_superbubble_annotations(Path(args.final_table), Path(args.open_superbubbles))
    metadata["openSuperbubbles"] = build_open_superbubble_lines(superbubbles)
    metadata["bubbles"] = build_bubble_annotations(Path(args.bubbles))
    metadata["render"] = {
        "raySteps": int(args.ray_steps),
        "opacityScale": float(args.opacity_scale),
        "brightness": float(args.brightness),
    }

    output_bin = Path(args.output_bin) if args.output_bin else None
    write_viewer(texture_order, metadata, Path(args.output_html), output_bin)
    print(f"interactive viewer HTML: {args.output_html}")
    if output_bin is not None:
        print(f"binary volume payload: {output_bin}")
    print(
        "volume texture: "
        f"{metadata['grid']['nx']} x {metadata['grid']['ny']} x {metadata['grid']['nz']} "
        f"= {metadata['grid']['totalVoxels']:,} cells"
    )
    print(
        "annotations: "
        f"{sum(len(g.get('items', [])) for g in metadata['openSuperbubbles'])} open superbubbles, "
        f"{len(metadata['bubbles'])} Grade A/B bubbles"
    )

def _format_runtime(seconds: float) -> str:
    """Return a compact human-readable wall-clock runtime string."""
    seconds = max(0.0, float(seconds))
    hours, remainder = divmod(seconds, 3600.0)
    minutes, seconds = divmod(remainder, 60.0)
    if hours >= 1.0:
        return f"{int(hours)} h {int(minutes):02d} min {seconds:05.2f} s"
    if minutes >= 1.0:
        return f"{int(minutes)} min {seconds:05.2f} s"
    return f"{seconds:.2f} s"


def _run_with_timing(entrypoint, script_file: str) -> None:
    """Run a script entry point and print its total wall-clock runtime."""
    script_name = Path(script_file).name
    start = time.perf_counter()
    try:
        entrypoint()
    except SystemExit as exc:
        elapsed = time.perf_counter() - start
        if exc.code in (None, 0):
            print(f"[runtime] {script_name} completed in: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        else:
            print(f"[runtime] {script_name} elapsed time before failure: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    except BaseException:
        elapsed = time.perf_counter() - start
        print(f"[runtime] {script_name} elapsed time before failure: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    elapsed = time.perf_counter() - start
    print(f"[runtime] {script_name} completed in: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)


if __name__ == "__main__":
    _run_with_timing(main, __file__)
