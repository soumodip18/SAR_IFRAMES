# SAR_IFRAMES
SAR_IFRAMES is the centralized repository for modular, embeddable dashboards and APIs that power the SAR Ventures ecosystem. It hosts self-contained IFrames for MRV (Measurement, Reporting & Verification), GIS mapping, aquaponics & hydroponics analytics, carbon accounting, and other regenerative agriculture insights. Each IFrame is designed to be plug-and-play, allowing internal teams, partners, and external stakeholders to embed live, data-driven dashboards directly into websites or portals.
Key Features:
Modular IFrames: Pre-built dashboards for various operational and analytical use-cases.
Centralized APIs: All data feeds powering the dashboards are versioned and maintained in one place.
Easy Integration: Embedding guide with query parameter support for farm ID, theme, and metrics.
Scalable & Secure: Version control, token-based access for sensitive MRV or carbon data, cloud-ready deployment.
Strategic Value: Provides transparency, auditability, and digital proof of SAR Ventures’ regenerative agriculture initiatives for partners, investors, and climate funds.
Intended Users:
Internal SAR Ventures teams (Mul Biotech, SequestraBionix, Digital Farm Solutions)
Partner organizations, NGOs, government agencies
Research collaborators and institutional investors
Impact:
This repository serves as the single source of truth for digital dashboards, turning field data into actionable insights, strengthening SAR Ventures’ ecosystem credibility, and creating a modular, enterprise-ready technology layer that supports expansion and monetization.
### Market Prices IFrame Dashboard

**Path:** `/dashboards/market_prices.html`  

**Description:** Dynamic IFrame dashboard showing variety-wise daily market prices. Users can select **market, commodity, and variety**, and the chart updates automatically.  

**Embed Example:**

```html
<iframe src="https://sarventures.github.io/iframes/market_prices.html" width="1000" height="700"></iframe>