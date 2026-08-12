"use client";

import { useEffect, useMemo, useRef, useState } from "react";

type Point = { x: number; y: number };
type ShapeName = "fillet" | "rotor" | "bracket";

type BallStep = {
  center: Point;
  radius: number;
  feature: Point;
};

const TAU = Math.PI * 2;

const shapeCopy: Record<ShapeName, { label: string; note: string }> = {
  fillet: {
    label: "Filleted block",
    note: "A superellipse stands in for a softened, machined CAD profile.",
  },
  rotor: {
    label: "Rotor profile",
    note: "A six-lobed section exposes repeated local feature sizes.",
  },
  bracket: {
    label: "Bracket section",
    note: "An asymmetric section creates a wider distribution of radii.",
  },
};

function length(point: Point) {
  return Math.hypot(point.x, point.y);
}

function normalise(point: Point): Point {
  const magnitude = length(point) || 1;
  return { x: point.x / magnitude, y: point.y / magnitude };
}

function pointForShape(shape: ShapeName, angle: number): Point {
  const cosine = Math.cos(angle);
  const sine = Math.sin(angle);

  if (shape === "fillet") {
    const exponent = 2 / 5.5;
    return {
      x: Math.sign(cosine) * Math.pow(Math.abs(cosine), exponent) * 1.02,
      y: Math.sign(sine) * Math.pow(Math.abs(sine), exponent) * 0.72,
    };
  }

  if (shape === "rotor") {
    const radius = 0.78 + 0.1 * Math.cos(6 * angle) + 0.035 * Math.cos(12 * angle);
    return { x: radius * cosine, y: radius * sine };
  }

  const radius = 0.72 + 0.13 * Math.cos(2 * angle) - 0.08 * Math.sin(3 * angle);
  return {
    x: radius * cosine + 0.13 * Math.cos(angle) ** 2,
    y: radius * sine * 0.9,
  };
}

function buildBoundary(shape: ShapeName, count: number): Point[] {
  return Array.from({ length: count }, (_, index) =>
    pointForShape(shape, (index / count) * TAU),
  );
}

function inwardNormal(points: Point[], index: number): Point {
  const previous = points[(index - 1 + points.length) % points.length];
  const next = points[(index + 1) % points.length];
  const tangent = normalise({ x: next.x - previous.x, y: next.y - previous.y });
  let normal = { x: -tangent.y, y: tangent.x };
  const point = points[index];

  if (normal.x * -point.x + normal.y * -point.y < 0) {
    normal = { x: -normal.x, y: -normal.y };
  }

  return normal;
}

function radiusForTangentPoint(points: Point[], index: number) {
  const point = points[index];
  const normal = inwardNormal(points, index);
  let radius = Number.POSITIVE_INFINITY;
  let feature = points[(index + Math.floor(points.length / 2)) % points.length];

  for (let candidateIndex = 0; candidateIndex < points.length; candidateIndex += 1) {
    const ringDistance = Math.min(
      Math.abs(candidateIndex - index),
      points.length - Math.abs(candidateIndex - index),
    );
    if (ringDistance < 3) continue;

    const candidate = points[candidateIndex];
    const delta = { x: candidate.x - point.x, y: candidate.y - point.y };
    const denominator = 2 * (normal.x * delta.x + normal.y * delta.y);
    if (denominator <= 0.0001) continue;

    const nextRadius = (delta.x * delta.x + delta.y * delta.y) / denominator;
    if (nextRadius > 0.015 && nextRadius < radius) {
      radius = nextRadius;
      feature = candidate;
    }
  }

  if (!Number.isFinite(radius)) radius = 0;
  return {
    radius,
    center: { x: point.x + normal.x * radius, y: point.y + normal.y * radius },
    feature,
    normal,
  };
}

function simulateShrinkingBall(
  points: Point[],
  index: number,
  initialRadius: number,
): BallStep[] {
  const point = points[index];
  const result = radiusForTangentPoint(points, index);
  const steps: BallStep[] = [];
  let radius = initialRadius;

  for (let iteration = 0; iteration < 8; iteration += 1) {
    const easedRadius = result.radius + (radius - result.radius) * 0.5;
    radius = Math.max(result.radius, easedRadius);
    steps.push({
      center: {
        x: point.x + result.normal.x * radius,
        y: point.y + result.normal.y * radius,
      },
      radius,
      feature: result.feature,
    });
    if (Math.abs(radius - result.radius) < 0.006) break;
  }

  steps.push({ center: result.center, radius: result.radius, feature: result.feature });
  return steps;
}

function formatNumber(value: number) {
  return value.toFixed(3);
}

function GeometryCanvas({
  points,
  candidates,
  selectedIndex,
  steps,
  threshold,
}: {
  points: Point[];
  candidates: ReturnType<typeof radiusForTangentPoint>[];
  selectedIndex: number;
  steps: BallStep[];
  threshold: number;
}) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const context = canvas.getContext("2d");
    if (!context) return;

    const render = () => {
      const bounds = canvas.getBoundingClientRect();
      const scaleFactor = window.devicePixelRatio || 1;
      canvas.width = Math.max(1, Math.round(bounds.width * scaleFactor));
      canvas.height = Math.max(1, Math.round(bounds.height * scaleFactor));
      context.setTransform(scaleFactor, 0, 0, scaleFactor, 0, 0);

      const width = bounds.width;
      const height = bounds.height;
      const scale = Math.min(width, height) * 0.35;
      const origin = { x: width * 0.5, y: height * 0.52 };
      const project = (point: Point): Point => ({
        x: origin.x + point.x * scale,
        y: origin.y - point.y * scale,
      });

      context.clearRect(0, 0, width, height);
      context.fillStyle = "#071713";
      context.fillRect(0, 0, width, height);

      context.strokeStyle = "rgba(211, 238, 227, 0.07)";
      context.lineWidth = 1;
      for (let x = origin.x % 32; x < width; x += 32) {
        context.beginPath();
        context.moveTo(x, 0);
        context.lineTo(x, height);
        context.stroke();
      }
      for (let y = origin.y % 32; y < height; y += 32) {
        context.beginPath();
        context.moveTo(0, y);
        context.lineTo(width, y);
        context.stroke();
      }

      context.beginPath();
      points.forEach((point, index) => {
        const projected = project(point);
        if (index === 0) context.moveTo(projected.x, projected.y);
        else context.lineTo(projected.x, projected.y);
      });
      context.closePath();
      context.fillStyle = "rgba(117, 240, 187, 0.055)";
      context.fill();
      context.strokeStyle = "#d8eee6";
      context.lineWidth = 1.4;
      context.stroke();

      candidates.forEach((candidate, index) => {
        if (index % 2 !== 0 || candidate.radius < threshold) return;
        const center = project(candidate.center);
        const strength = Math.min(1, candidate.radius / 0.7);
        context.fillStyle = `rgba(117, 240, 187, ${0.25 + strength * 0.55})`;
        context.beginPath();
        context.arc(center.x, center.y, 1.25 + strength * 1.6, 0, TAU);
        context.fill();
      });

      const point = project(points[selectedIndex]);
      const current = steps[steps.length - 1];
      steps.slice(0, -1).forEach((step, index) => {
        const center = project(step.center);
        context.strokeStyle = `rgba(255, 177, 95, ${0.08 + index * 0.025})`;
        context.lineWidth = 1;
        context.beginPath();
        context.arc(center.x, center.y, step.radius * scale, 0, TAU);
        context.stroke();
      });

      if (current) {
        const center = project(current.center);
        const feature = project(current.feature);
        context.fillStyle = "rgba(255, 177, 95, 0.08)";
        context.strokeStyle = "#ffb15f";
        context.lineWidth = 1.6;
        context.beginPath();
        context.arc(center.x, center.y, current.radius * scale, 0, TAU);
        context.fill();
        context.stroke();

        context.setLineDash([5, 5]);
        context.strokeStyle = "rgba(255, 177, 95, 0.68)";
        context.beginPath();
        context.moveTo(point.x, point.y);
        context.lineTo(feature.x, feature.y);
        context.stroke();
        context.setLineDash([]);

        context.fillStyle = "#ffb15f";
        context.beginPath();
        context.arc(center.x, center.y, 4, 0, TAU);
        context.fill();
        context.fillStyle = "#f5fff9";
        context.beginPath();
        context.arc(point.x, point.y, 4.5, 0, TAU);
        context.fill();
        context.fillStyle = "#ffb15f";
        context.beginPath();
        context.arc(feature.x, feature.y, 3.5, 0, TAU);
        context.fill();
      }
    };

    render();
    const resizeObserver = new ResizeObserver(render);
    resizeObserver.observe(canvas);
    return () => resizeObserver.disconnect();
  }, [candidates, points, selectedIndex, steps, threshold]);

  return (
    <canvas
      ref={canvasRef}
      className="geometry-canvas"
      aria-label="Animated two-dimensional reconstruction of the medial-axis shrinking-ball algorithm"
      role="img"
    />
  );
}

export default function Home() {
  const [shape, setShape] = useState<ShapeName>("fillet");
  const [detail, setDetail] = useState(160);
  const [threshold, setThreshold] = useState(0.18);
  const [initialRadius, setInitialRadius] = useState(1.05);
  const [phase, setPhase] = useState(0.12);
  const [playing, setPlaying] = useState(true);

  const points = useMemo(() => buildBoundary(shape, detail), [detail, shape]);
  const candidates = useMemo(
    () => points.map((_, index) => radiusForTangentPoint(points, index)),
    [points],
  );
  const selectedIndex = Math.floor(phase * points.length) % points.length;
  const steps = useMemo(
    () => simulateShrinkingBall(points, selectedIndex, initialRadius),
    [initialRadius, points, selectedIndex],
  );
  const { accepted, medianRadius } = useMemo(() => {
    const retained = candidates.filter((candidate) => candidate.radius >= threshold);
    const sortedRadii = retained.map((candidate) => candidate.radius).sort((a, b) => a - b);
    return {
      accepted: retained,
      medianRadius: sortedRadii.length ? sortedRadii[Math.floor(sortedRadii.length / 2)] : 0,
    };
  }, [candidates, threshold]);
  const selectedRadius = candidates[selectedIndex]?.radius ?? 0;

  const histogram = useMemo(() => {
    const bins = Array.from({ length: 14 }, () => 0);
    candidates.forEach((candidate) => {
      const bin = Math.min(bins.length - 1, Math.floor(candidate.radius * bins.length));
      bins[Math.max(0, bin)] += 1;
    });
    return bins;
  }, [candidates]);
  const maxBin = Math.max(...histogram, 1);

  useEffect(() => {
    if (!playing || window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;
    let frame = 0;
    let previous = performance.now();
    const animate = (time: number) => {
      if (time - previous > 55) {
        setPhase((value) => (value + 0.0028) % 1);
        previous = time;
      }
      frame = requestAnimationFrame(animate);
    };
    frame = requestAnimationFrame(animate);
    return () => cancelAnimationFrame(frame);
  }, [playing]);

  return (
    <main>
      <nav className="topbar" aria-label="Project navigation">
        <a className="monogram" href="#top" aria-label="Back to project introduction">
          CG<span>/21</span>
        </a>
        <div className="nav-links">
          <a href="#method">Method</a>
          <a href="#contribution">Contribution</a>
          <a href="https://github.com/FLABDUL/beng-thesis" target="_blank" rel="noreferrer">
            Source
          </a>
        </div>
      </nav>

      <section className="hero" id="top">
        <div className="hero-copy">
          <p className="eyebrow">BEng thesis / computational geometry / 2021</p>
          <h1>Development of a machine-learning CAD filter using computational geometry.</h1>
          <p className="lede">
            An experimental C++ pipeline that samples triangulated CAD surfaces,
            approximates their medial axis, and exposes local feature radii for
            downstream filtering. This browser reconstruction makes the geometry tangible.
          </p>
          <div className="tech-line" aria-label="Project technologies">
            <span>C++</span><span>PCL</span><span>Eigen</span><span>NumPy</span><span>Jupyter</span>
          </div>
        </div>
        <aside className="hero-aside" aria-label="Project summary">
          <p className="aside-number">01</p>
          <p>
            Start with an oriented point cloud. Roll a virtual sphere inward from each
            surface point until a second feature constrains it.
          </p>
          <dl>
            <div><dt>Input</dt><dd>CAD mesh</dd></div>
            <div><dt>Signal</dt><dd>Medial radius</dd></div>
            <div><dt>Output</dt><dd>Filter features</dd></div>
          </dl>
        </aside>
      </section>

      <section className="lab" aria-labelledby="lab-title">
        <div className="section-heading">
          <p className="eyebrow">Interactive reconstruction</p>
          <h2 id="lab-title">Roll the shrinking ball.</h2>
          <p>
            The orange sphere begins at the selected surface point. Its centre travels
            along the inward normal until another feature point limits the radius.
          </p>
        </div>

        <div className="lab-grid">
          <div className="viewport">
            <div className="viewport-meta">
              <span>Section / {shapeCopy[shape].label}</span>
              <span className={selectedRadius >= threshold ? "accepted" : "rejected"}>
                {selectedRadius >= threshold ? "Feature retained" : "Feature filtered"}
              </span>
            </div>
            <GeometryCanvas
              points={points}
              candidates={candidates}
              selectedIndex={selectedIndex}
              steps={steps}
              threshold={threshold}
            />
            <div className="legend" aria-label="Diagram legend">
              <span><i className="dot dot-surface" />surface</span>
              <span><i className="dot dot-axis" />retained medial centres</span>
              <span><i className="dot dot-ball" />active ball</span>
            </div>
          </div>

          <div className="controls">
            <fieldset>
              <legend>Profile</legend>
              <div className="segmented">
                {(Object.keys(shapeCopy) as ShapeName[]).map((shapeName) => (
                  <button
                    type="button"
                    className={shape === shapeName ? "active" : ""}
                    aria-pressed={shape === shapeName}
                    onClick={() => setShape(shapeName)}
                    key={shapeName}
                  >
                    {shapeCopy[shapeName].label}
                  </button>
                ))}
              </div>
              <p className="control-note">{shapeCopy[shape].note}</p>
            </fieldset>

            <div className="range-control">
              <span><b id="initial-radius-label">Initial ball radius</b><output>{initialRadius.toFixed(2)}</output></span>
              <input
                id="initial-radius"
                aria-labelledby="initial-radius-label"
                type="range" min="0.72" max="1.35" step="0.01"
                value={initialRadius}
                onChange={(event) => setInitialRadius(Number(event.target.value))}
              />
            </div>

            <div className="range-control">
              <span><b id="filter-threshold-label">Filter threshold</b><output>{threshold.toFixed(2)}</output></span>
              <input
                id="filter-threshold"
                aria-labelledby="filter-threshold-label"
                type="range" min="0.04" max="0.55" step="0.01"
                value={threshold}
                onChange={(event) => setThreshold(Number(event.target.value))}
              />
            </div>

            <div className="range-control">
              <span><b id="surface-samples-label">Surface samples</b><output>{detail}</output></span>
              <input
                id="surface-samples"
                aria-labelledby="surface-samples-label"
                type="range" min="96" max="256" step="8"
                value={detail}
                onChange={(event) => setDetail(Number(event.target.value))}
              />
            </div>

            <div className="transport">
              <button type="button" onClick={() => setPlaying((value) => !value)}>
                {playing ? "Pause scan" : "Resume scan"}
              </button>
              <button
                type="button"
                className="secondary"
                onClick={() => {
                  setPlaying(false);
                  setPhase((value) => (value + 0.035) % 1);
                }}
              >
                Step point
              </button>
            </div>

            <div className="readout" aria-live="polite">
              <div><span>Active radius</span><strong>{formatNumber(selectedRadius)}</strong></div>
              <div><span>Iterations</span><strong>{steps.length - 1}</strong></div>
              <div><span>Feature index</span><strong>{selectedIndex.toString().padStart(3, "0")}</strong></div>
            </div>
          </div>
        </div>

        <div className="metrics">
          <article><p>Input samples</p><strong>{points.length}</strong><span>oriented boundary points</span></article>
          <article><p>Retained centres</p><strong>{accepted.length}</strong><span>above radius threshold</span></article>
          <article><p>Median radius</p><strong>{formatNumber(medianRadius)}</strong><span>normalised section units</span></article>
          <article className="histogram-card">
            <div className="histogram-heading"><p>Radius distribution</p><span>0 to 1</span></div>
            <div className="histogram" aria-label="Histogram of medial radii">
              {histogram.map((count, index) => (
                <i key={index} style={{ height: `${Math.max(7, (count / maxBin) * 100)}%` }} />
              ))}
            </div>
          </article>
        </div>
      </section>

      <section className="method" id="method" aria-labelledby="method-title">
        <div className="section-heading sticky-heading">
          <p className="eyebrow">The method</p>
          <h2 id="method-title">From solid model to geometric signal.</h2>
        </div>
        <ol className="process-list">
          <li><span>01</span><div><h3>Sample the surface</h3><p>Subdivide the STL mesh and export vertices with their outward normals as contiguous NumPy arrays.</p></div></li>
          <li><span>02</span><div><h3>Shrink tangent balls</h3><p>For every oriented point, search the point cloud and iteratively move a sphere centre along its normal.</p></div></li>
          <li><span>03</span><div><h3>Recover local radii</h3><p>Record the converged medial centre, constraining feature index, and ball radius introduced by this fork.</p></div></li>
          <li><span>04</span><div><h3>Filter CAD features</h3><p>Use the radius distribution as a scale-aware signal for separating fine geometric details from the primary form.</p></div></li>
        </ol>
      </section>

      <section className="contribution" id="contribution" aria-labelledby="contribution-title">
        <div>
          <p className="eyebrow">What I contributed</p>
          <h2 id="contribution-title">Research code, made observable.</h2>
        </div>
        <div className="contribution-copy">
          <p>
            The project began with the open-source <a href="https://github.com/tudelft3d/masbcpp" target="_blank" rel="noreferrer">masbcpp</a> implementation.
            My thesis fork adapted its I/O and shrinking-ball path for a CAD experiment,
            increased convergence precision, and surfaced the radius of every computed ball.
          </p>
          <ul>
            <li>Added inner and outer medial-radius NumPy outputs.</li>
            <li>Built an STL to remesh to normals to MAT notebook workflow.</li>
            <li>Instrumented convergence and feature-point behaviour for analysis.</li>
            <li>Explored radius histograms as a CAD feature-filtering signal.</li>
          </ul>
          <p className="provenance">
            The visual demo is a two-dimensional educational reconstruction. The production
            research implementation remains the C++/PCL pipeline in the repository.
          </p>
        </div>
      </section>

      <footer>
        <p>Abdul Hakim Norazman / BEng project revitalised for the web</p>
        <a href="https://github.com/FLABDUL/beng-thesis" target="_blank" rel="noreferrer">Explore the C++ source</a>
      </footer>
    </main>
  );
}
