class LanguageSwitcher {
    constructor(config = {}) {
        this.languages = config.languages || [
            { code: 'en', name: 'English', default: true },
            { code: 'zh_CN', name: '中文' }
        ];
        this.currentLang = this.detectCurrentLanguage();
        this.init();
    }

    // Detect current language from URL path
    detectCurrentLanguage() {
        const path = window.location.pathname;
        const langCodes = this.languages.map(lang => lang.code);
        
        for (const lang of langCodes) {
            if (path.includes(`/${lang}/`)) {
                return lang;
            }
        }
        
        // Fallback to default language
        const defaultLang = this.languages.find(lang => lang.default) || this.languages[0];
        return defaultLang.code;
    }

    // Get relative path without language prefix
    getRelativePath() {
        const currentPath = window.location.pathname;
        const currentLang = this.currentLang;
        
        // Remove the language segment from path
        const langSegment = `/${currentLang}/`;
        const langIndex = currentPath.indexOf(langSegment);
        
        if (langIndex !== -1) {
            return currentPath.substring(langIndex + langSegment.length - 1); // Keep leading slash
        }
        
        // If no language segment found, return current path
        return currentPath;
    }

    // Generate URL for a specific language
    generateLangUrl(targetLang) {
        const currentUrl = new URL(window.location.href);
        const relativePath = this.getRelativePath();
        
        // Replace or add language segment
        const basePath = currentUrl.pathname.split('/').slice(0, -1).join('/');
        const newPath = `/${targetLang}${relativePath}`;
        
        currentUrl.pathname = newPath;
        return currentUrl.toString();
    }

    // Create switcher HTML
    createSwitcherHTML() {
        const container = document.createElement('div');
        container.className = 'sphinx-language-switcher';
        
        const select = document.createElement('select');
        select.className = 'language-select';
        select.innerHTML = this.languages.map(lang => 
            `<option value="${lang.code}" ${this.currentLang === lang.code ? 'selected' : ''}>
                ${lang.name}
            </option>`
        ).join('');

        select.addEventListener('change', (e) => {
            const targetLang = e.target.value;
            window.location.href = this.generateLangUrl(targetLang);
        });

        container.appendChild(select);
        return container;
    }

    // Auto-inject into page
    init() {
        // Wait for DOM to be ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', () => this.inject());
        } else {
            this.inject();
        }
    }

    inject() {
        // Try to find common Sphinx header locations
        const headerSelectors = [
            '.related:first-child',
            '.sphinxsidebar .sidebar-container:first-child',
            '.body .section:first-child',
            'header',
            '.navbar'
        ];

        for (const selector of headerSelectors) {
            const element = document.querySelector(selector);
            if (element) {
                const switcher = this.createSwitcherHTML();
                element.prepend(switcher);
                break;
            }
        }

        // Fallback: insert at top of body
        if (!document.querySelector('.sphinx-language-switcher')) {
            const switcher = this.createSwitcherHTML();
            document.body.insertBefore(switcher, document.body.firstChild);
        }
    }
}

// Auto-initialize
document.addEventListener('DOMContentLoaded', function() {
    new LanguageSwitcher();
});